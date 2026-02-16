#version 330 core //#version 300 es
precision highp float;
uvec4 R_STATE;

in vec2 vTexCoord;
out vec4 fragColor;

// main uniforms
uniform sampler2D u_sample;
uniform float u_sample_part;

uniform float u_time;
uniform vec2 u_resolution;
uniform sampler2D u_skybox;
uniform vec2 u_mouse;
uniform vec3 u_pos;

uniform vec2 u_seed1;
uniform vec2 u_seed2;
uniform bool RT; // false == lite mode

// object data arrays
#define MAX_LEN 20
uniform float len;
uniform float type[MAX_LEN];
uniform vec4 color[MAX_LEN];
uniform float mat[MAX_LEN];
uniform vec3 pos[MAX_LEN];
uniform vec3 size[MAX_LEN];
uniform float extra[MAX_LEN];

const int MAX_REF = 8;

struct rayHit{
    vec2 dist;
    vec4 color;
    vec3 normal;
    int mat;
};

vec2 sphIntersect(in vec3 ro, in vec3 rd, in vec3 ce, float ra);
float groundIntersect(in vec3 ro, in vec3 rd, in float h);
float IntersectCylinder(vec3 ro, vec3 rd, vec3 a, vec3 b, float r);
vec3 computeCylinderNormal(vec3 a, vec3 b, vec3 pos);
vec2 boxIntersection(in vec3 ro, in vec3 rd, in vec3 boxPos, in vec3 boxSize, inout vec3 outNormal);

void materialRay(inout vec3 ro, inout vec3 rd, in rayHit ray);
rayHit RayAttributes(in vec3 ro, in vec3 rd, in rayHit ray, in vec3 light, in int ind);
rayHit ObjectAttributes(in vec3 ro, in vec3 rd, in rayHit ray, in int ind);
rayHit castRay(in vec3 ro, in vec3 rd, in vec3 light, out int ind);

uint TausStep(uint z, int S1, int S2, int S3, uint M);
uint LCGStep(uint z, uint A, uint C);
float random();
vec3 randomOnSphere();
vec2 hash22(vec2 p);

// for mouse controls
mat2 rot(float a) {
	float s = sin(a);
	float c = cos(a);
	return mat2(c, -s, s, c);
}

vec3 skybox(vec3 rd, in vec3 light){
    return vec3(0.);
    if(rd.z > 0.9) rd.z = 0.9 - (rd.z - 0.9);
    vec2 uv = vec2(atan(rd.x, abs(rd.y)), -asin(abs(rd.z))*4.+1.9);
    uv /= 3.14159265;
    uv = uv* 0.5 + 0.5;
    vec3 col = (texture(u_skybox, uv).rgb);
    if (length(col) > length(vec3(0.7))) col *= (length(col) - length(vec3(0.7))) + 1.;
    vec3 sun = vec3(1.0, 1.0, 0.8);
    sun *= max(0.0, pow(dot(rd, light), 2048.0))*1000.;
    return sun + col*0.8;
}

vec3 rayShader(in vec3 ro, in vec3 rd){
    vec3 light = normalize(vec3(-0.5, -0.5, 0.5));
    int ind=0;
    rayHit inter;
    vec3 color = vec3(1.);
    for(int i=0; i < MAX_REF; i++){
        inter = castRay(ro, rd, light, ind);
        inter = RayAttributes(ro, rd, inter, light, ind);
        inter = ObjectAttributes(ro, rd, inter, ind);
        materialRay(ro, rd, inter);
        color *= abs(inter.color.xyz);
        if(rd == vec3(0.)){break;};
    }
    return color;
}

vec3 ContrastSaturationBrightness(vec3 color, float brt, float sat, float con)
{
	// increase or decrease theese values to adjust r, g and b color channels seperately
	const float AvgLumR = 0.5;
	const float AvgLumG = 0.5;
	const float AvgLumB = 0.5;
	
	const vec3 LumCoeff = vec3(0.2125, 0.7154, 0.0721);
	
	vec3 AvgLumin  = vec3(AvgLumR, AvgLumG, AvgLumB);
	vec3 brtColor  = color * brt;
	vec3 intensity = vec3(dot(brtColor, LumCoeff));
	vec3 satColor  = mix(intensity, brtColor, sat);
	vec3 conColor  = mix(AvgLumin, satColor, con);
	
	return conColor;
}

vec3 colorThreshold(vec3 color) {
    float maxValue = max(max(color.r, color.g), color.b);
    if (maxValue > 1.0) {
        color /= maxValue; 
    }
    return color;
}

void main() {
    // converting texture coordinates from [0,1] into [-1,1]
    vec2 uv = vTexCoord * 2.0 - 2.;
    vec2 uvRes = hash22(uv + 1.0) * u_resolution + u_resolution;
    // normalizing the ratio
    uv.x *= u_resolution.x / u_resolution.y;

    vec3 ro = u_pos;

    // shooting rays from the center of the screen using uv
    // FOV can be corrected here
    vec3 rd = normalize(vec3(1.0, uv / 4.));
    rd.xz *= rot(u_mouse.y);
	rd.xy *= rot(-u_mouse.x);
    
    // random function seeds
    R_STATE.x = uint(u_seed1.x + uvRes.x);
	R_STATE.y = uint(u_seed1.y + uvRes.x);
	R_STATE.z = uint(u_seed2.x + uvRes.y);
	R_STATE.w = uint(u_seed2.y + uvRes.y);

    // render cycle
	int samples = RT? 10 : 1; // samples
    vec3 col = vec3(0.0);
    for(int i = 0; i < samples; i++) {
        col += rayShader(ro, rd);
    }
    col /= float(samples);

    // optional color correction
    /* float exposure = 1.0;
    float white = 1.0;
	col *= white * exposure;
	col = (col * (1.0 + col / white / white)) / (1.0 + col);
    col = ContrastSaturationBrightness(col, 1., 1., 1.0); */

    col = colorThreshold(col);

    vec3 sampleCol = texture(u_sample, vTexCoord/2.0).rgb;
    col = mix(sampleCol, col, u_sample_part);
    fragColor = vec4(col.xyz, 1.);
}

// materials. modifies ray direction and origin
void materialRay(inout vec3 ro, inout vec3 rd, in rayHit ray){
    switch(ray.mat){
        default:{ // no light calculations. use for light sources
            ro = ro + rd * ray.dist.x + ray.normal * 0.00001;
            rd = vec3(0.);
            break;
        }
        case 1:{ // mirror
            ro = ro + rd * ray.dist.x + ray.normal * 0.00001;
            rd = normalize(reflect(rd, ray.normal));
            break;
        }
        case 2:{ // glass
            if(dot(rd, ray.normal) > -1. + fract(random())*0.5){
                ro = ro + rd * ray.dist.x + 0.00001;
                rd = normalize(reflect(rd, ray.normal));
            } else{
                ro = ro + rd * ray.dist.y + 0.00001;
                rd = refract(rd, ray.normal, (1.0 / (1.0 - ray.color.a))); // refraction factor
            }
            //vec3 refracted = refract(rd, ray.normal, (1.0 / (1.0 - ray.color.a))); // refraction factor
            //vec3 reflected = normalize(reflect(rd, ray.normal));
            break;
        }
        case 3:{ // matte with roughness
            ro = ro + rd * ray.dist.x + ray.normal * 0.00001;
            vec3 reflected = reflect(rd, ray.normal);
            vec3 r = randomOnSphere();
	        vec3 diffuse = normalize(r * dot(r, ray.normal));
	        rd = mix(diffuse, reflected, ray.color.a);
            break;
        }
    }
}

// custom object attributes based on index
rayHit ObjectAttributes(in vec3 ro, in vec3 rd, in rayHit ray, in int ind){
    switch(ind){
        case 4:{
            vec3 pos = vec3(ray.dist.x*rd + ro)*0.5+2.;
            float diff = abs(pos.z - tan(pos.x)-4.5);
            ray.color.xyz = vec3(0.5 + 0.5 * sin(diff), 0.5 + 0.5 * sin(diff + 2.0), 0.5 + 0.5 * sin(diff + 4.0))*2.5+1.25;
            break;
        }
        case 7:{
            vec3 pos = vec3(ray.dist.x*rd + ro)*3.;
            pos.x += 0.3;
            ray.color.xyz = vec3(0.5 + 0.5 * sin(pos.x), 0.5 + 0.5 * sin(pos.x + 2.), 0.5 + 0.5 * sin(pos.x + 4.))*128.5+1.25;
            break;
        }
        case 8:{
            vec3 pos = vec3(ray.dist.x*rd + ro)*0.75;
            float diff = abs(pos.z - tan(pos.y)-4.5);
            ray.color.xyz = vec3(0.5 + 0.5 * sin(diff), 0.5 + 0.5 * sin(diff + 2.), 0.5 + 0.5 * sin(diff + 4.))*5.5+1.25;
            break;
        }
    }
    return ray;
}

// object type attributes
rayHit RayAttributes(in vec3 ro, in vec3 rd, in rayHit ray, in vec3 light, in int ind){
    vec3 outNormal = vec3(1.);
    switch(ind >= 0? int(type[ind]) : ind){
        case -2: { // sky
            ray.normal = vec3(0.0);
            ray.color.xyz = skybox(rd, light);
            ray.color.a = 1.;
            //rd = vec3(0.);
            return ray;
            break;
        }
        case -1: { // ground
            vec3 grid = vec3(ray.dist.x*rd + ro);
            if((floor(mod(grid.x+0.125, 4.)/0.2) == 0.) || (floor(mod(grid.y+0.125, 4.)/0.2) == 0.)){
                ray.color = vec4(0.4, 0.4, 0.7, 0.);
            } else {ray.color = vec4(0.7, 0.7, 0.7, 0.);}
            ray.normal = vec3(0., 0., 1.);
            break;
        }
        case 1: { // sphere
            ray.normal = normalize(ro + rd*ray.dist.x - pos[ind]);
            break;
        }
        case 2: { // box
            ray.dist = vec2(boxIntersection(ro, rd, pos[ind], size[ind], outNormal)); // VERY INEFFICIENT
            ray.normal = outNormal;
            break;
        }
        case 3: { // cylinder from 2 points and a radius
            ray.normal = computeCylinderNormal(pos[ind], size[ind], ro + rd * ray.dist.x + ray.normal * 0.00001);
            break;
        }
    }
    return ray;
}

// intersection check on all objects in the scene
rayHit castRay(in vec3 ro, in vec3 rd, in vec3 light, out int ind){
    rayHit ray;

    // search for closest intersected object
    vec2 intersect;
    intersect.xy = vec2(1. / 0.);
    vec3 outNormal=vec3(0.);
    ind=-2;
    for(int i=0; i<int(len); i++){
        vec2 new_inter;
        switch(int(type[i])){
            case 1: new_inter=sphIntersect(ro, rd, pos[i], size[i].x); break;
            case 2: new_inter=vec2(boxIntersection(ro, rd, pos[i], size[i], outNormal)); break;
            case 3: new_inter=vec2(IntersectCylinder(ro, rd, pos[i], size[i], extra[i])); break;
            default: new_inter=vec2(1. / 0.);
        }
        if(new_inter.x < intersect.x && new_inter.x > 0.){intersect=new_inter; ind = i;}
    }

    ray.mat = int(mat[ind]);

    // sky and ground
    vec2 new_inter = vec2(groundIntersect(ro, rd, 0.));
    if(new_inter.x < intersect.x && new_inter.x > 0.){intersect=new_inter; ind = -1; ray.mat = RT? 3 : 0;}

    ray.dist = intersect;
    ray.color = color[ind];

    return ray;
}

// object intersection functions

vec2 sphIntersect(in vec3 ro, in vec3 rd, in vec3 ce, float ra) {
    vec3 oc = ro - ce;
    float b = dot(oc, rd);
    float c = dot(oc, oc) - ra*ra;
    float h = b*b - c;
    if (h < 0.0) return vec2(-1.0);
    h = sqrt(h);
    return vec2(-b - h, -b + h);
}

float groundIntersect(in vec3 ro, in vec3 rd, in float h){
    if(rd.z > 0. || ro.z < h) return -1.;
    vec3 it = (((abs(ro.z-h))/rd.z)*rd);
    return length(it);
}

// don't even ask me what is up with cylinder
float IntersectCylinder(vec3 ro, vec3 rd, vec3 a, vec3 b, float r) {
    vec3 ba = b - a;
    vec3 oa = ro - a;
    float baba = dot(ba, ba);
    float bard = dot(ba, rd);
    float baoa = dot(ba, oa);

    vec3 rdPerp = rd * baba - ba * bard;
    vec3 oaPerp = oa * baba - ba * baoa;
    float a_ = dot(rdPerp, rdPerp);
    float b_ = dot(rdPerp, oaPerp);
    float c_ = dot(oaPerp, oaPerp) - r*r * baba*baba;
    float h = b_*b_ - a_*c_;

    float t = 1e20;

    // Боковая поверхность
    if (h > 0.0) {
        float t0 = (-b_ - sqrt(h)) / a_;
        if (t0 > 0.0) {
            float y = baoa + t0 * bard;
            if (y > 0.0 && y < baba)
                t = t0;
        }
    }

    // Крышки
    vec3 n = normalize(ba);
    float denom = dot(rd, n);
    if (abs(denom) > 1e-6) {
        // Нижняя крышка
        float t1 = dot(a - ro, n) / denom;
        if (t1 > 0.0 && length(ro + rd * t1 - a) <= r && t1 < t)
            t = t1;
        // Верхняя крышка
        float t2 = dot(b - ro, n) / denom;
        if (t2 > 0.0 && length(ro + rd * t2 - b) <= r && t2 < t)
            t = t2;
    }

    return (t < 1e19) ? t-1e-5 : -1.0;
}

vec3 computeCylinderNormal(vec3 a, vec3 b, vec3 pos) {
     // Вектор направления цилиндра
    vec3 cylinderAxis = normalize(b - a);
    
    // Вектор от основания цилиндра до точки pos
    vec3 toPos = pos - a;
    
    // Проекция toPos на ось цилиндра
    float projectionLength = dot(toPos, cylinderAxis);
    vec3 projection = projectionLength * cylinderAxis;
    
    // Вектор от проекции до точки pos
    vec3 radialVector = toPos - projection;

    // Проверка, находится ли pos на боковой поверхности или на крышках
    if (projectionLength < 0.0) {
        // Нижняя крышка
        return normalize(vec3(0.0, 0.0, -1.0)); // Нормаль направлена вниз
    } else if (projectionLength > length(b - a)) {
        // Верхняя крышка
        return normalize(vec3(0.0, 0.0, 1.0)); // Нормаль направлена вверх
    } else {
        // Боковая поверхность
        return normalize(radialVector); // Нормаль направлена от оси цилиндра
    }
}

vec2 boxIntersection( in vec3 ro, in vec3 rd, in vec3 boxPos, in vec3 boxSize, inout vec3 outNormal ) 
{
    vec3 localRo = ro - boxPos;

    vec3 m = 1.0 / rd;
    vec3 n = m * localRo;
    vec3 k = abs(m) * boxSize;
    vec3 t1 = -n - k;
    vec3 t2 = -n + k;
    
    float tN = max(max(t1.x, t1.y), t1.z);
    float tF = min(min(t2.x, t2.y), t2.z);

    if( tN > tF || tF < 0.0 ) return vec2(-1.0);
    if(outNormal.x!=0.){
        // Determine the normal by checking which component of t1 equals tN
        if( tN > 0.0 ) {
            if( tN == t1.x ) outNormal = vec3(-sign(rd.x), 0.0, 0.0);
            else if( tN == t1.y ) outNormal = vec3(0.0, -sign(rd.y), 0.0);
            else outNormal = vec3(0.0, 0.0, -sign(rd.z));
        } else {
            if( tF == t2.x ) outNormal = vec3(sign(rd.x), 0.0, 0.0);
            else if( tF == t2.y ) outNormal = vec3(0.0, sign(rd.y), 0.0);
            else outNormal = vec3(0.0, 0.0, sign(rd.z));
        }
    }
    return vec2(tN, tF);
}

// functions for random

uint TausStep(uint z, int S1, int S2, int S3, uint M)
{
	uint b = (((z << S1) ^ z) >> S2);
	return (((z & M) << S3) ^ b);	
}

uint LCGStep(uint z, uint A, uint C)
{
	return (A * z + C);	
}

float random()
{
	R_STATE.x = TausStep(R_STATE.x, 13, 19, 12, uint(49467294));
	R_STATE.y = TausStep(R_STATE.y, 2, 25, 4, uint(42497288));
	R_STATE.z = TausStep(R_STATE.z, 3, 11, 17, uint(42967280));
	R_STATE.w = LCGStep(R_STATE.w, uint(1664525), uint(1013904223));
	return 2.3283064365387e-10 * float((R_STATE.x ^ R_STATE.y ^ R_STATE.z ^ R_STATE.w));
}

vec3 randomOnSphere() {
	vec3 rand = vec3(random(), random(), random());
	float theta = rand.x * 2.0 * 3.14159265;
	float v = rand.y;
	float phi = acos(2.0 * v - 1.0);
	float r = pow(rand.z, 1.0 / 3.0);
	float x = r * sin(phi) * cos(theta);
	float y = r * sin(phi) * sin(theta);
	float z = r * cos(phi);
	return vec3(x, y, z);
}

vec2 hash22(vec2 p)
{
	p += u_seed1.x;
	vec3 p3 = fract(vec3(p.xyx) * vec3(.1031, .1030, .0973));
	p3 += dot(p3, p3.yzx+33.33);
	return fract((p3.xx+p3.yz)*p3.zy);
}