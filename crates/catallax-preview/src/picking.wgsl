// Unproved nearest named triangle in a field chart preview.
struct Triangle { a: vec4<f32>, b: vec4<f32>, c: vec4<f32> }
@group(0) @binding(0) var<storage,read> points: array<vec4<f32>>;
@group(0) @binding(1) var<storage,read> triangles: array<Triangle>;
@group(0) @binding(2) var<storage,read_write> selected: array<vec4<u32>>;
fn segment_distance(p: vec3<f32>, a: vec3<f32>, b: vec3<f32>) -> f32 {
    let edge=b-a;
    let t=clamp(dot(p-a,edge)/dot(edge,edge),0.0,1.0);
    let d=p-(a+t*edge);
    return dot(d,d);
}
fn triangle_distance(p: vec3<f32>, t: Triangle) -> f32 {
    let a=t.a.xyz; let b=t.b.xyz; let c=t.c.xyz;
    let normal=cross(b-a,c-a);
    let n=normal/sqrt(dot(normal,normal));
    let value=dot(p-a,n);
    let foot=p-n*(value);
    if(dot(cross(b-a,foot-a),n)>=0.0 &&
       dot(cross(c-b,foot-b),n)>=0.0 &&
       dot(cross(a-c,foot-c),n)>=0.0) {
        return value*value;
    }
    return min(segment_distance(p,a,b),min(segment_distance(p,b,c),segment_distance(p,c,a)));
}
@compute @workgroup_size(64)
fn select_triangles(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x>=arrayLength(&points)) {return;}
    var best=bitcast<f32>(0x7f800000u); var owner=0u;
    for(var i=0u;i<arrayLength(&triangles);i+=1u) {
        let d=triangle_distance(points[id.x].xyz,triangles[i]);
        if(d<best) {best=d;owner=i+1u;}
    }
    selected[id.x]=vec4<u32>(owner,bitcast<u32>(best),0u,0u);
}
