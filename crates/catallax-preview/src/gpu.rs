//! The graphics-card device shared by the unproved preview routes.
use catallax_numbers::Refusal;
use wgpu::util::DeviceExt;

/// A wgpu adapter and device for unproved preview work. It owns no proof
/// state and its results are display values only.
pub struct Gpu {
    pub(crate) device: wgpu::Device,
    pub(crate) queue: wgpu::Queue,
    name: String,
}

impl Gpu {
    /// Opens the high-performance adapter. `disabled` selects no device, so
    /// callers can exercise their processor route deliberately. The error is
    /// the reason code: `preview.gpu-disabled`, `preview.no-adapter: ...` or
    /// `preview.no-device: ...`.
    pub fn new(disabled: bool) -> Result<Self, String> {
        if disabled {
            return Err("preview.gpu-disabled".into());
        }
        pollster::block_on(async {
            let instance = wgpu::Instance::new(&wgpu::InstanceDescriptor::default());
            let adapter = instance
                .request_adapter(&wgpu::RequestAdapterOptions {
                    power_preference: wgpu::PowerPreference::HighPerformance,
                    ..Default::default()
                })
                .await
                .map_err(|e| format!("preview.no-adapter: {e}"))?;
            let info = adapter.get_info();
            let (device, queue) = adapter
                .request_device(&wgpu::DeviceDescriptor::default())
                .await
                .map_err(|e| format!("preview.no-device: {e}"))?;
            Ok(Self {
                device,
                queue,
                name: format!("{} {:?} {:?}", info.name, info.backend, info.device_type),
            })
        })
    }

    /// Adapter name, backend and device type.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Compiles a WGSL compute pipeline. A validation failure is the typed
    /// refusal `preview.shader`, never a silent fallback.
    pub(crate) fn pipeline(
        &self,
        source: &str,
        entry: &str,
    ) -> Result<wgpu::ComputePipeline, Refusal> {
        self.device.push_error_scope(wgpu::ErrorFilter::Validation);
        let shader = self
            .device
            .create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("unproved field preview"),
                source: wgpu::ShaderSource::Wgsl(source.into()),
            });
        let pipeline = self
            .device
            .create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
                label: Some("unproved preview"),
                layout: None,
                module: &shader,
                entry_point: Some(entry),
                compilation_options: Default::default(),
                cache: None,
            });
        if let Some(error) = pollster::block_on(self.device.pop_error_scope()) {
            return Err(Refusal::new("preview.shader", error.to_string()));
        }
        Ok(pipeline)
    }

    /// Waits for submitted work and maps `readback`. Device loss, a timeout or
    /// a failed mapping is the typed refusal `preview.gpu-device`.
    pub(crate) fn read<T: bytemuck::Pod>(
        &self,
        readback: &wgpu::Buffer,
        bytes: u64,
    ) -> Result<Vec<T>, Refusal> {
        let device = |e: String| Refusal::new("preview.gpu-device", e);
        let slice = readback.slice(..bytes);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |r| {
            let _ = tx.send(r);
        });
        self.device
            .poll(wgpu::PollType::Wait {
                submission_index: None,
                timeout: Some(std::time::Duration::from_secs(120)),
            })
            .map_err(|e| device(e.to_string()))?;
        rx.recv()
            .map_err(|e| device(e.to_string()))?
            .map_err(|e| device(e.to_string()))?;
        let mapped = slice.get_mapped_range();
        let result = bytemuck::cast_slice(&mapped).to_vec();
        drop(mapped);
        readback.unmap();
        Ok(result)
    }

    /// Runs one entry point once over fresh storage buffers bound in order,
    /// and reads back buffer `output`. Used by the throwaway viewer.
    pub fn dispatch<T: bytemuck::Pod>(
        &self,
        source: &str,
        entry: &str,
        data: &[&[u8]],
        output: usize,
        groups: [u32; 3],
    ) -> Result<Vec<T>, Refusal> {
        let pipeline = self.pipeline(source, entry)?;
        let buffers: Vec<_> = data
            .iter()
            .map(|bytes| {
                self.device
                    .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                        label: Some("unproved preview buffer"),
                        contents: bytes,
                        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
                    })
            })
            .collect();
        let buffer = &buffers[output];
        let readback = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("unproved preview readback"),
            size: buffer.size(),
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        let group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &pipeline.get_bind_group_layout(0),
            entries: &buffers
                .iter()
                .enumerate()
                .map(|(i, buffer)| wgpu::BindGroupEntry {
                    binding: i as u32,
                    resource: buffer.as_entire_binding(),
                })
                .collect::<Vec<_>>(),
        });
        let mut encoder = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = encoder.begin_compute_pass(&Default::default());
            pass.set_pipeline(&pipeline);
            pass.set_bind_group(0, &group, &[]);
            pass.dispatch_workgroups(groups[0], groups[1], groups[2]);
        }
        encoder.copy_buffer_to_buffer(buffer, 0, &readback, 0, buffer.size());
        self.queue.submit([encoder.finish()]);
        self.read(&readback, buffer.size())
    }
}
