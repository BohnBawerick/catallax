// Unproved display only. Each invocation owns one pixel.
@group(0) @binding(0) var<storage, read_write> pixels: array<vec4<f32>>;
const WIDTH: u32 = 960u;
const HEIGHT: u32 = 720u;
const SCALE: f32 = 0.065;
const TARGET: vec3<f32> = vec3<f32>(20.0,15.0,10.0);
const FORWARD: vec3<f32> = vec3<f32>(-0.577350269,-0.577350269,-0.577350269);
const RIGHT: vec3<f32> = vec3<f32>(0.707106781,-0.707106781,0.0);
const UP: vec3<f32> = vec3<f32>(-0.408248290,-0.408248290,0.816496581);
@compute @workgroup_size(64)
fn probe(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x>=arrayLength(&pixels)) { return; }
    let p=pixels[id.x].xyz;
    pixels[id.x]=vec4<f32>(p,field(p));
}
@compute @workgroup_size(8,8)
fn render(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x>=WIDTH || id.y>=HEIGHT) { return; }
    let index=id.y*WIDTH+id.x;
    let origin=TARGET-FORWARD*80.0+RIGHT*(f32(id.x)+0.5-f32(WIDTH)*0.5)*SCALE+UP*(f32(HEIGHT)*0.5-f32(id.y)-0.5)*SCALE;
    let a=(vec3<f32>(-0.1)-origin)/FORWARD;
    let b=(vec3<f32>(40.1,30.1,22.1)-origin)/FORWARD;
    let lo=min(a,b); let hi=max(a,b);
    var t=max(0.0,max(lo.x,max(lo.y,lo.z)));
    let end=min(hi.x,min(hi.y,hi.z));
    var previous=t;
    pixels[index]=vec4<f32>(0.0);
    for(var step=0u;step<16384u;step+=1u) {
        if(t>end) { return; }
        let p=origin+t*FORWARD;
        let d=field(p);
        if(d<0.0) {
            var l=previous; var h=t;
            for(var refine=0u;refine<16u;refine+=1u) {
                let m=(l+h)*0.5;
                if(field(origin+m*FORWARD)<0.0) { h=m; } else { l=m; }
            }
            pixels[index]=vec4<f32>(origin+((l+h)*0.5)*FORWARD,1.0);
            return;
        }
        previous=t;
        t+=max(d*0.8,0.002);
    }
    pixels[index]=vec4<f32>(0.0,0.0,0.0,-1.0);
}
