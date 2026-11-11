// Unproved nearest named closed patch of the fixed P3 recipe.
@group(0) @binding(0) var<storage, read> hits: array<vec4<f32>>;
@group(0) @binding(1) var<storage, read_write> face_ids: array<u32>;
fn square(x: f32) -> f32 { return x*x; }
fn rect_distance(p: vec3<f32>, lo: vec3<f32>, hi: vec3<f32>) -> f32 {
    let q=p-clamp(p,lo,hi); return dot(q,q);
}
fn rim_at(p: vec3<f32>, extent: f32, angle: f32) -> f32 {
    let d=1.0-cos(angle);
    let q=vec3<f32>(d,clamp(p.y,d,extent-d),9.0+sin(angle));
    return dot(p-q,p-q);
}
fn rim_distance(p: vec3<f32>, extent: f32) -> f32 {
    let half_pi=1.5707963267948966;
    var best=1e30; var index=0u;
    for(var i=0u;i<=32u;i+=1u) {
        let d=rim_at(p,extent,half_pi*f32(i)/32.0);
        if(d<best) {best=d;index=i;}
    }
    var lo=half_pi*f32(max(i32(index)-1,0))/32.0;
    var hi=half_pi*f32(min(index+1u,32u))/32.0;
    for(var i=0u;i<28u;i+=1u) {
        let a=(2.0*lo+hi)/3.0; let b=(lo+2.0*hi)/3.0;
        if(rim_at(p,extent,a)<rim_at(p,extent,b)) {hi=b;} else {lo=a;}
    }
    return min(best,rim_at(p,extent,(lo+hi)*0.5));
}
fn owner(p: vec3<f32>) -> u32 {
    let rho=length(p.xy-vec2<f32>(20.0,15.0));
    let blend_radius=8.0+JUNCTION_RADIUS;
    var top=clamp(p.xy,vec2<f32>(1.0),vec2<f32>(39.0,29.0));
    if(length(top-vec2<f32>(20.0,15.0))<blend_radius) {
        var direction=vec2<f32>(1.0,0.0);
        if(rho>0.0) {direction=(p.xy-vec2<f32>(20.0,15.0))/rho;}
        top=vec2<f32>(20.0,15.0)+blend_radius*direction;
    }
    let delta=p.xy-top;
    let torus=vec2<f32>(rho-blend_radius,p.z-(10.0+JUNCTION_RADIUS));
    var junction=min(dot(torus+vec2<f32>(JUNCTION_RADIUS,0.0),torus+vec2<f32>(JUNCTION_RADIUS,0.0)),dot(torus+vec2<f32>(0.0,JUNCTION_RADIUS),torus+vec2<f32>(0.0,JUNCTION_RADIUS)));
    if(torus.x<=0.0 && torus.y<=0.0) {junction=square(length(torus)-JUNCTION_RADIUS);}
    let distances=array<f32,13>(
        rect_distance(p,vec3<f32>(0.0),vec3<f32>(40.0,30.0,0.0)),
        dot(delta,delta)+square(p.z-10.0),
        rect_distance(p,vec3<f32>(0.0),vec3<f32>(0.0,30.0,9.0)),
        rect_distance(p,vec3<f32>(40.0,0.0,0.0),vec3<f32>(40.0,30.0,9.0)),
        rect_distance(p,vec3<f32>(0.0),vec3<f32>(40.0,0.0,9.0)),
        rect_distance(p,vec3<f32>(0.0,30.0,0.0),vec3<f32>(40.0,30.0,9.0)),
        square(max(rho-8.0,0.0))+square(p.z-22.0),
        square(rho-8.0)+square(p.z-clamp(p.z,10.0+JUNCTION_RADIUS,22.0)),
        junction,
        rim_distance(vec3<f32>(p.x,p.y,p.z),30.0),
        rim_distance(vec3<f32>(40.0-p.x,p.y,p.z),30.0),
        rim_distance(vec3<f32>(p.y,p.x,p.z),40.0),
        rim_distance(vec3<f32>(30.0-p.y,p.x,p.z),40.0)
    );
    var best=1e30; var result=0u;
    for(var i=0u;i<13u;i+=1u) {
        if(distances[i]<best) {best=distances[i];result=NAMED_FACES[i];}
    }
    return result;
}
@compute @workgroup_size(64)
fn select_faces(@builtin(global_invocation_id) id: vec3<u32>) {
    if(id.x>=arrayLength(&hits)) {return;}
    face_ids[id.x]=0u;
    if(hits[id.x].w==1.0) {face_ids[id.x]=owner(hits[id.x].xyz);}
}
