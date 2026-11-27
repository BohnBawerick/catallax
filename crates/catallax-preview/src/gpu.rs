//! Compute on the renderer's device and queue. The preview never maps a buffer.

use wgpu::util::DeviceExt;

use crate::mesher::{GpuRefusal, MeshBuffers, PreviewVertex};

const SHADER: &str = concat!(
    include_str!("mesher/compute.wgsl"),
    include_str!("mesher/analytic.wgsl"),
    include_str!("mesher/grid.wgsl"),
);
const GROUP_SIZE: u32 = 64;

pub(crate) struct Gpu {
    device: wgpu::Device,
    queue: wgpu::Queue,
    pipeline: wgpu::ComputePipeline,
    layout: wgpu::BindGroupLayout,
    pub(crate) storage_bytes: usize,
    pub(crate) index_bytes: usize,
    groups_per_axis: u32,
}

impl Gpu {
    pub(crate) async fn new(
        device: &wgpu::Device,
        queue: &wgpu::Queue,
    ) -> Result<Self, GpuRefusal> {
        let limits = device.limits();
        if limits.max_storage_buffers_per_shader_stage < 4
            || limits.max_compute_invocations_per_workgroup < GROUP_SIZE
            || limits.max_compute_workgroup_size_x < GROUP_SIZE
            || limits.max_compute_workgroups_per_dimension == 0
        {
            return Err(GpuRefusal::UnsupportedDevice);
        }
        let scope = device.push_error_scope(wgpu::ErrorFilter::Validation);
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("Catallax preview mesher"),
            source: wgpu::ShaderSource::Wgsl(SHADER.into()),
        });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("Catallax preview mesher"),
            layout: None,
            module: &shader,
            entry_point: Some("mesh"),
            compilation_options: Default::default(),
            cache: None,
        });
        if let Some(error) = scope.pop().await {
            return Err(GpuRefusal::Shader(error.to_string()));
        }
        let layout = pipeline.get_bind_group_layout(0);
        Ok(Self {
            device: device.clone(),
            queue: queue.clone(),
            pipeline,
            layout,
            storage_bytes: limits
                .max_buffer_size
                .min(limits.max_storage_buffer_binding_size)
                .min(u64::from(u32::MAX)) as usize,
            index_bytes: limits.max_buffer_size.min(u32::MAX as u64) as usize,
            groups_per_axis: limits.max_compute_workgroups_per_dimension,
        })
    }

    pub(crate) fn max_vertices(&self) -> usize {
        (self.storage_bytes / size_of::<PreviewVertex>()).min(
            (u64::from(self.groups_per_axis).pow(2) * u64::from(GROUP_SIZE)).min(u32::MAX as u64)
                as usize,
        )
    }

    pub(crate) async fn dispatch(
        &self,
        inputs: &[u8],
        positions: &[u8],
        surfaces: &[u8],
        indices: &[u8],
        vertex_count: u32,
    ) -> Result<MeshBuffers, GpuRefusal> {
        let oom = self.device.push_error_scope(wgpu::ErrorFilter::OutOfMemory);
        let validation = self.device.push_error_scope(wgpu::ErrorFilter::Validation);
        let input = self.upload("Catallax preview UVs", inputs, wgpu::BufferUsages::STORAGE);
        let positions = self.upload(
            "Catallax preview edge and grid positions",
            positions,
            wgpu::BufferUsages::STORAGE,
        );
        let surfaces = self.upload(
            "Catallax preview surface parameters",
            surfaces,
            wgpu::BufferUsages::STORAGE,
        );
        let index = self.upload(
            "Catallax preview indices",
            indices,
            wgpu::BufferUsages::INDEX | wgpu::BufferUsages::COPY_SRC,
        );
        let vertex = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Catallax preview vertices"),
            size: u64::from(vertex_count) * size_of::<PreviewVertex>() as u64,
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::VERTEX
                | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });
        let buffers = [&input, &positions, &surfaces, &vertex];
        let group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &self.layout,
            entries: &std::array::from_fn::<_, 4, _>(|binding| wgpu::BindGroupEntry {
                binding: binding as u32,
                resource: buffers[binding].as_entire_binding(),
            }),
        });
        let mut encoder = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = encoder.begin_compute_pass(&Default::default());
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &group, &[]);
            let groups = vertex_count.div_ceil(GROUP_SIZE);
            let x = groups.min(self.groups_per_axis);
            pass.dispatch_workgroups(x, groups.div_ceil(x), 1);
        }
        self.queue.submit([encoder.finish()]);
        let validation = validation.pop();
        let oom = oom.pop();
        if let Some(error) = validation.await.or(oom.await) {
            return Err(GpuRefusal::Device(error.to_string()));
        }
        Ok(MeshBuffers {
            vertex,
            index,
            vertex_count,
        })
    }

    fn upload(&self, label: &str, contents: &[u8], usage: wgpu::BufferUsages) -> wgpu::Buffer {
        self.device
            .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: Some(label),
                contents,
                usage,
            })
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn mesher_shaders_validate() {
        let module = naga::front::wgsl::parse_str(super::SHADER)
            .unwrap_or_else(|error| panic!("{}", error.emit_to_string(super::SHADER)));
        naga::valid::Validator::new(
            naga::valid::ValidationFlags::all(),
            naga::valid::Capabilities::empty(),
        )
        .validate(&module)
        .unwrap();
    }
}
