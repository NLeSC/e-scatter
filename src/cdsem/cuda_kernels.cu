/**
 * @file src/cdsem/cuda_kernels.cu
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include "cuda_kernels.cuh"
#include <cfloat>
#include <common/cuda_make_ptr.cuh>
#include <common/cuda_vec3_math.cuh>

__device__ const float eps = 10.0f*FLT_EPSILON;
__device__ const float mc2 = 5.1099897e+05f;

/*
namespace fermi_sea {
__device__ const int n = 32;
__device__ const float z1 = 1.0000000e-03f;
__device__ const float z2 = 1.0000000e+01f;
__device__ const float p1 = 2.1088175e-05f;
__device__ const float p2 = 5.4595325e+01f;
__device__ const float z2p[n] = {
    -1.0766798e+01f, -1.0321033e+01f, -9.8752317e+00f, -9.4293823e+00f,
    -8.9834681e+00f, -8.5374670e+00f, -8.0913477e+00f, -7.6450720e+00f,
    -7.1985850e+00f, -6.7518148e+00f, -6.3046651e+00f, -5.8570080e+00f,
    -5.4086747e+00f, -4.9594412e+00f, -4.5090156e+00f, -4.0570192e+00f,
    -3.6029665e+00f, -3.1462471e+00f, -2.6861079e+00f, -2.2216477e+00f,
    -1.7518265e+00f, -1.2755001e+00f, -7.9148793e-01f, -2.9867163e-01f,
    +2.0388594e-01f, +7.1682256e-01f, +1.2403889e+00f, +1.7744150e+00f,
    +2.3183398e+00f, +2.8712943e+00f, +3.4322128e+00f, +3.9999483e+00f
};
__device__ const float p2z[n] = {
    -6.9077554e+00f, -6.5902653e+00f, -6.2728038e+00f, -5.9553800e+00f,
    -5.6380086e+00f, -5.3207092e+00f, -5.0035086e+00f, -4.6864424e+00f,
    -4.3695607e+00f, -4.0529299e+00f, -3.7366409e+00f, -3.4208145e+00f,
    -3.1056123e+00f, -2.7912462e+00f, -2.4779902e+00f, -2.1661916e+00f,
    -1.8562782e+00f, -1.5487578e+00f, -1.2442060e+00f, -9.4323510e-01f,
    -6.4644718e-01f, -3.5437316e-01f, -6.7411400e-02f, +2.1421888e-01f,
    +4.9049795e-01f, +7.6159561e-01f, +1.0278367e+00f, +1.2896539e+00f,
    +1.5475388e+00f, +1.8019993e+00f, +2.0535290e+00f, +2.3025851e+00f
};
};
*/

template<typename T>
__device__ T clamp(T x, T x1, T x2) {
    return max(x1, min(x2, x));
}

__device__ float interp1(const float* ptr, int pitch, int height, int ix, int iy, int iz, float sx) {
    ptr = cuda_make_ptr(ptr, pitch, height, iy, iz);
    return (1.0f-sx)*ptr[ix]+sx*ptr[ix+1];
}

__device__ float interp2(const float* ptr, int pitch, int height, int ix, int iy, int iz, float sx, float sy) {
    return sy*interp1(ptr, pitch, height, ix, iy, iz, sx)+(1.0f-sy)*interp1(ptr, pitch, height, ix, iy+1, iz, sx);
}

__global__ void __init_rand_state(curandState* rand_state_p, unsigned long long seed, int n) {
    const int i = threadIdx.x+blockIdx.x*blockDim.x;
    if(i >= n)
        return;
    curand_init(seed, i, 0, &(rand_state_p[i]));
}

__global__ void __init_trajectory(cuda_particle_struct pstruct, cuda_geometry_struct gstruct, cuda_material_struct mstruct, curandState* rand_state_dev_p) {
    const int pid = threadIdx.x+blockIdx.x*blockDim.x;
    if(pid >= pstruct.capacity)
        return;

    int status = pstruct.status_dev_p[pid];
    if(status == cuda_particle_struct::TERMINATED)
        return;
    else if(status == cuda_particle_struct::DETECTED)
        return;
    else if(status == cuda_particle_struct::PENDING)
        return;

    const auto __inside_octree = [&](const float3& pos, const float3& center, const float3& size) {
        if((pos.x > center.x-0.5f*size.x) && (pos.x < center.x+0.5f*size.x))
        if((pos.y > center.y-0.5f*size.y) && (pos.y < center.y+0.5f*size.y))
        if((pos.z > center.z-0.5f*size.z) && (pos.z < center.z+0.5f*size.z))
            return true;
        return false;
    };

    const float3 pos = make_float3(
        pstruct.pos_x_dev_p[pid],
        pstruct.pos_y_dev_p[pid],
        pstruct.pos_z_dev_p[pid]
    );

    status = cuda_particle_struct::TERMINATED;

    if(__inside_octree(pos, gstruct.root_center, gstruct.root_size)) {
        const int mid = pstruct.material_idx_dev_p[pid];
        float distance;
        if(mid < 0) {
            // vacuum propgation
            status = cuda_particle_struct::NO_EVENT;
            distance = norm3df(gstruct.root_size.x, gstruct.root_size.y, gstruct.root_size.z);
        } else {
            // material propagation
            if(pstruct.K_energy_dev_p[pid] < mstruct.barrier_dev_p[mid]) {
                // EXIT: if energy is below barrier
                status = cuda_particle_struct::TERMINATED;
            } else {
                // determine attenuation length
                curandState rand_state = rand_state_dev_p[pid];

                const int& pitch = mstruct.table_pitch;
                const float& K1 = mstruct.K_energy_range.x;
                const float& K2 = mstruct.K_energy_range.y;
                const int& nK = mstruct.table_dim.x;
                const int& nP = mstruct.table_dim.y;
                const float x = __fdividef(__logf(pstruct.K_energy_dev_p[pid]/K1), __logf(K2/K1))*(nK-1);
                const int ix = clamp(__float2int_rd(x), 0, nK-2);

                const float elastic_imfp = __expf(interp1(mstruct.elastic_dev_p, pitch, nP+1, ix, 0, mid, x-ix));
                const float inelastic_imfp = __expf(interp1(mstruct.inelastic_dev_p, pitch, nP+1, ix, 0, mid, x-ix));
                const float total_imfp = elastic_imfp+inelastic_imfp;
                distance = __fdividef(-__logf(curand_uniform(&rand_state)), total_imfp);
                status = cuda_particle_struct::ELASTIC_EVENT;
                if(curand_uniform(&rand_state) > elastic_imfp/total_imfp)
                    status = cuda_particle_struct::INELASTIC_EVENT;
                rand_state_dev_p[pid] = rand_state;
            }
        }
        pstruct.distance_dev_p[pid] = distance;
    }
    pstruct.status_dev_p[pid] = status;
}

__global__ void __update_trajectory(cuda_particle_struct pstruct, cuda_geometry_struct gstruct, cuda_material_struct mstruct) {
    const int pid = threadIdx.x+blockIdx.x*blockDim.x;
    if(pid >= pstruct.capacity)
        return;

    int status = pstruct.status_dev_p[pid];
    if(status == cuda_particle_struct::TERMINATED)
        return;
    else if(status == cuda_particle_struct::DETECTED)
        return;
    else if(status == cuda_particle_struct::PENDING)
        return;

    float3 pos = make_float3(
        pstruct.pos_x_dev_p[pid],
        pstruct.pos_y_dev_p[pid],
        pstruct.pos_z_dev_p[pid]
    );
    float3 dir = make_float3(
        pstruct.dir_x_dev_p[pid],
        pstruct.dir_y_dev_p[pid],
        pstruct.dir_z_dev_p[pid]
    );
    dir = dir*rnorm3df(dir.x, dir.y, dir.z);
    const float K = pstruct.K_energy_dev_p[pid];
    const int mid = pstruct.material_idx_dev_p[pid];
    const int tid = pstruct.triangle_idx_dev_p[pid];

    float distance = pstruct.distance_dev_p[pid];
    uint64_t location = 1;
    do {
        // traverse to location
        float3 center = gstruct.root_center;
        float3 size = gstruct.root_size;
        int index = 0;
        for(int i = 20-__clzll(location)/3; i >= 0; i--) {
            const int octant = (location>>(i*3))&7;
            index = cuda_make_ptr(gstruct.octree_dev_p, gstruct.octree_pitch, index)[octant];
            center.x += 0.25f*size.x*(2.0f*(octant&1)-1.0f);
            center.y += 0.25f*size.y*(1.0f*(octant&2)-1.0f);
            center.z += 0.25f*size.z*(0.5f*(octant&4)-1.0f);
            size.x *= 0.5f;
            size.y *= 0.5f;
            size.z *= 0.5f;
        }

        // traverse to leaf
        while(index >= 0) {
            int octant = 0;
            octant += (pos.x > center.x) ? 1 : 0;
            octant += (pos.y > center.y) ? 2 : 0;
            octant += (pos.z > center.z) ? 4 : 0;
            location = (location<<3)|octant;
            index = cuda_make_ptr(gstruct.octree_dev_p, gstruct.octree_pitch, index)[octant];
            center.x += 0.25f*size.x*(2.0f*(octant&1)-1.0f);
            center.y += 0.25f*size.y*(1.0f*(octant&2)-1.0f);
            center.z += 0.25f*size.z*(0.5f*(octant&4)-1.0f);
            size.x *= 0.5f;
            size.y *= 0.5f;
            size.z *= 0.5f;
        }

        // determine intersections
        float intersect;
        int target_id;
        const float tx = __fdividef(center.x+copysignf(0.5f*size.x+eps, dir.x)-pos.x, dir.x);
        const float ty = __fdividef(center.y+copysignf(0.5f*size.y+eps, dir.y)-pos.y, dir.y);
        const float tz = __fdividef(center.z+copysignf(0.5f*size.z+eps, dir.z)-pos.z, dir.z);
        if((tx < ty) && (tx < tz)) {
            intersect = tx;
            target_id = -1;
        } else if((ty < tx) && (ty < tz)) {
            intersect = ty;
            target_id = -2;
        } else {
            intersect = tz;
            target_id = -4;
        }
        for(int i = 0; i < gstruct.occupancy; i++) {
            const int j = cuda_make_ptr(gstruct.octree_dev_p, gstruct.octree_pitch, -index)[i];
            if(j < 0)
                break;
            if(j == tid)
                continue;
            const float3 T = make_float3(
                gstruct.triangle_Ax_dev_p[j],
                gstruct.triangle_Ay_dev_p[j],
                gstruct.triangle_Az_dev_p[j]
            );
            const float3 e1 = make_float3(
                gstruct.triangle_Bx_dev_p[j]-T.x,
                gstruct.triangle_By_dev_p[j]-T.y,
                gstruct.triangle_Bz_dev_p[j]-T.z
            );
            const float3 e2 = make_float3(
                gstruct.triangle_Cx_dev_p[j]-T.x,
                gstruct.triangle_Cy_dev_p[j]-T.y,
                gstruct.triangle_Cz_dev_p[j]-T.z
            );
            int mat_idx_in = gstruct.material_idx_in_dev_p[j];
            int mat_idx_out = gstruct.material_idx_out_dev_p[j];
            if(dot_product(cross_product(e1, e2), dir) < 0)
                mat_idx_out = mat_idx_in;
            if((mat_idx_out == mid) || (mat_idx_out == cuda_geometry_struct::NOP))
                continue;
            else if((mat_idx_out == cuda_geometry_struct::DETECTOR_LT50) && (K >= 50))
                continue;
            else if((mat_idx_out == cuda_geometry_struct::DETECTOR_GE50) && (K < 50))
                continue;
            // T. Möller and B. Trumbore, Journal of Graphics Tools, 2(1):21--28, 1997.
            const float3 pvec = cross_product(dir, e2);
            const float det = dot_product(e1, pvec);
            if((det > -eps) && (det < eps))
                continue;
            const float3 tvec = pos-T;
            const float u = __fdividef(dot_product(tvec, pvec), det);
            if((u < -eps) || (u > 1.0f+eps))
                continue;
            const float3 qvec = cross_product(tvec, e1);
            const float v = __fdividef(dot_product(dir, qvec), det);
            if((v < -eps) || (u+v > 1.0f+eps))
                continue;
            const float t = __fdividef(dot_product(e2, qvec), det);
            if((t > 0) && (t <= intersect+eps)) {
                intersect = t;
                target_id = j;
            }
        }

        // manage intersections
        if(intersect >= distance) {
            // EXIT: no intersection
            pstruct.pos_x_dev_p[pid] = pos.x+dir.x*distance;
            pstruct.pos_y_dev_p[pid] = pos.y+dir.y*distance;
            pstruct.pos_z_dev_p[pid] = pos.z+dir.z*distance;
            return;
        } else if(target_id >= 0) {
            // EXIT: triangle intersection
            pstruct.status_dev_p[pid] = cuda_particle_struct::INTERSECT_EVENT;
            pstruct.triangle_idx_dev_p[pid] = target_id;
            pstruct.pos_x_dev_p[pid] = pos.x+dir.x*intersect;
            pstruct.pos_y_dev_p[pid] = pos.y+dir.y*intersect;
            pstruct.pos_z_dev_p[pid] = pos.z+dir.z*intersect;
            return;
        }
        distance -= intersect;
        pos.x += dir.x*intersect;
        pos.y += dir.y*intersect;
        pos.z += dir.z*intersect;

        // find adjacent node
        unsigned int mask = -target_id;
        unsigned int value;
        if(mask == 1)
            value = (dir.x >= 0) ? 0 : 1;
        else if(mask == 2)
            value = (dir.y >= 0) ? 0 : 2;
        else
            value = (dir.z >= 0) ? 0 : 4;
        while(location > 1) {
            if((location&mask) == value) {
                location ^= mask;
                break;
            }
            location >>= 3;
        }
    } while(location > 1);

    // EXIT: out of grid
    pstruct.status_dev_p[pid] = cuda_particle_struct::TERMINATED;
    return;
}

__global__ void __apply_isec_event(cuda_particle_struct pstruct, cuda_geometry_struct gstruct, cuda_material_struct mstruct, curandState* rand_state_dev_p) {
    const int i = threadIdx.x+blockIdx.x*blockDim.x;
    if(i >= pstruct.capacity)
        return;
    const int pid = pstruct.particle_idx_dev_p[i];

    int status = pstruct.status_dev_p[pid];
    if(status != cuda_particle_struct::INTERSECT_EVENT)
        return;

    float3 dir = make_float3(
        pstruct.dir_x_dev_p[pid],
        pstruct.dir_y_dev_p[pid],
        pstruct.dir_z_dev_p[pid]
    );
    dir = dir*rnorm3df(dir.x, dir.y, dir.z);

    const int tid = pstruct.triangle_idx_dev_p[pid];
    const float3 T = make_float3(
        gstruct.triangle_Ax_dev_p[tid],
        gstruct.triangle_Ay_dev_p[tid],
        gstruct.triangle_Az_dev_p[tid]
    );
    const float3 e1 = make_float3(
        gstruct.triangle_Bx_dev_p[tid]-T.x,
        gstruct.triangle_By_dev_p[tid]-T.y,
        gstruct.triangle_Bz_dev_p[tid]-T.z
    );
    const float3 e2 = make_float3(
        gstruct.triangle_Cx_dev_p[tid]-T.x,
        gstruct.triangle_Cy_dev_p[tid]-T.y,
        gstruct.triangle_Cz_dev_p[tid]-T.z
    );
    float3 normal = cross_product(e1, e2);
    normal = normal*rnorm3df(normal.x, normal.y, normal.z);
    const float cos_alpha = dot_product(normal, dir);

    int mid_in, mid_out;
    if(cos_alpha > 0) {
        mid_in = gstruct.material_idx_in_dev_p[tid];
        mid_out = gstruct.material_idx_out_dev_p[tid];
    } else {
        mid_in = gstruct.material_idx_out_dev_p[tid];
        mid_out = gstruct.material_idx_in_dev_p[tid];
    }
    switch(mid_out) {
        case cuda_geometry_struct::DETECTOR:
        case cuda_geometry_struct::DETECTOR_LT50:
        case cuda_geometry_struct::DETECTOR_GE50:
            pstruct.status_dev_p[pid] = cuda_particle_struct::DETECTED;
            return;
        case cuda_geometry_struct::TERMINATOR:
            pstruct.status_dev_p[pid] = cuda_particle_struct::TERMINATED;
            return;
        case cuda_geometry_struct::MIRROR:
            pstruct.dir_x_dev_p[pid] = dir.x-2.0f*normal.x*cos_alpha;
            pstruct.dir_y_dev_p[pid] = dir.y-2.0f*normal.y*cos_alpha;
            pstruct.dir_z_dev_p[pid] = dir.z-2.0f*normal.z*cos_alpha;
            return;
        default:
            break;
    }

    float dU = 0;
    if(mid_out >= 0)
        dU += mstruct.barrier_dev_p[mid_out];
    if(mid_in >= 0)
        dU -= mstruct.barrier_dev_p[mid_in];

    // R. Shimizu and Z. J. Ding, Rep. Prog. Phys., 55, 487-531, 1992
    //  see Eqs. 3.20, 3.23 and 3.24
    curandState rand_state = rand_state_dev_p[pid];
    const float K = pstruct.K_energy_dev_p[pid];
    const float z = sqrtf(1.0f+dU/(K*cos_alpha*cos_alpha));
    if((K*cos_alpha*cos_alpha+dU > 0) && (curand_uniform(&rand_state) < __fdividef(4.0f*z, ((1.0f+z)*(1.0f+z))))) {
        pstruct.dir_x_dev_p[pid] = (dir.x-normal.x*cos_alpha)+normal.x*cos_alpha*z;
        pstruct.dir_y_dev_p[pid] = (dir.y-normal.y*cos_alpha)+normal.y*cos_alpha*z;
        pstruct.dir_z_dev_p[pid] = (dir.z-normal.z*cos_alpha)+normal.z*cos_alpha*z;
        pstruct.K_energy_dev_p[pid] = K+dU;
        pstruct.material_idx_dev_p[pid] = mid_out;
    } else if((dU < 0) && (curand_uniform(&rand_state) < __expf(1.0f+0.5f*K/dU))) {
        // surface absorption? (see Kieft & Bosch code)
        pstruct.status_dev_p[pid] = cuda_particle_struct::TERMINATED;
    } else {
        // total internal reflection
        pstruct.dir_x_dev_p[pid] = dir.x-2.0f*normal.x*cos_alpha;
        pstruct.dir_y_dev_p[pid] = dir.y-2.0f*normal.y*cos_alpha;
        pstruct.dir_z_dev_p[pid] = dir.z-2.0f*normal.z*cos_alpha;
    }
    rand_state_dev_p[pid] = rand_state;
}

__global__ void __apply_elastic_event(cuda_particle_struct pstruct, cuda_material_struct mstruct, curandState* rand_state_dev_p) {
    const int i = threadIdx.x+blockIdx.x*blockDim.x;
    if(i >= pstruct.capacity)
        return;
    const int pid = pstruct.particle_idx_dev_p[i];

    int status = pstruct.status_dev_p[pid];
    if(status != cuda_particle_struct::ELASTIC_EVENT)
        return;

    pstruct.triangle_idx_dev_p[pid] = -1;

    curandState rand_state = rand_state_dev_p[pid];
    const float x = __fdividef(__logf(pstruct.K_energy_dev_p[pid]/mstruct.K_energy_range.x), __logf(mstruct.K_energy_range.y/mstruct.K_energy_range.x))*(mstruct.table_dim.x-1);
    const int ix = clamp(__float2int_rd(x), 0, mstruct.table_dim.x-2);
    const float y = curand_uniform(&rand_state)*(mstruct.table_dim.y-1);
    const int iy = clamp(__float2int_rd(y), 0, mstruct.table_dim.y-2);
    const int mid = pstruct.material_idx_dev_p[pid];
    const float cos_theta = clamp(interp2(mstruct.elastic_dev_p, mstruct.table_pitch, mstruct.table_dim.y+1, ix, 1+iy, mid, x-ix, y-iy), -1.0f, 1.0f);
    const float sin_theta = sqrtf(1.0f-cos_theta*cos_theta);
    float3 dir = make_float3(pstruct.dir_x_dev_p[pid], pstruct.dir_y_dev_p[pid], pstruct.dir_z_dev_p[pid]);
    dir = dir*rnorm3df(dir.x, dir.y, dir.z);
    float sin_azimuth, cos_azimuth;
    __sincosf(atan2f(dir.y, dir.x), &sin_azimuth, &cos_azimuth);
    const float3 unit_v = make_float3(dir.z*cos_azimuth, dir.z*sin_azimuth, -sqrtf(__saturatef(1.0f-dir.z*dir.z)));
    const float3 unit_u = cross_product(unit_v, dir);
    float sin_phi, cos_phi;
    sincospif(2.0f*curand_uniform(&rand_state), &sin_phi, &cos_phi);
    pstruct.dir_x_dev_p[pid] = dir.x*cos_theta+(unit_u.x*cos_phi+unit_v.x*sin_phi)*sin_theta;
    pstruct.dir_y_dev_p[pid] = dir.y*cos_theta+(unit_u.y*cos_phi+unit_v.y*sin_phi)*sin_theta;
    pstruct.dir_z_dev_p[pid] = dir.z*cos_theta+(unit_u.z*cos_phi+unit_v.z*sin_phi)*sin_theta;
    rand_state_dev_p[pid] = rand_state;
}

__global__ void __apply_inelastic_event(cuda_particle_struct pstruct, cuda_material_struct mstruct, curandState* rand_state_dev_p) {
    const int i = threadIdx.x+blockIdx.x*blockDim.x;
    if(i >= pstruct.capacity)
        return;
    const int pid = pstruct.particle_idx_dev_p[i];

    int status = pstruct.status_dev_p[pid];
    if(status != cuda_particle_struct::INELASTIC_EVENT)
    if(status != cuda_particle_struct::PENDING)
        return;

    const int sid = pstruct.particle_idx_dev_p[pstruct.capacity-1-i];
    if(pstruct.status_dev_p[sid] != cuda_particle_struct::TERMINATED) {
        pstruct.status_dev_p[pid] = cuda_particle_struct::PENDING;
        return;
    }
    if(status == cuda_particle_struct::PENDING)
        pstruct.status_dev_p[pid] = cuda_particle_struct::INELASTIC_EVENT;

    pstruct.triangle_idx_dev_p[pid] = -1;

    curandState rand_state = rand_state_dev_p[pid];
    const float K = pstruct.K_energy_dev_p[pid];
    const float x = __fdividef(__logf(pstruct.K_energy_dev_p[pid]/mstruct.K_energy_range.x), __logf(mstruct.K_energy_range.y/mstruct.K_energy_range.x))*(mstruct.table_dim.x-1);
    const int ix = clamp(__float2int_rd(x), 0, mstruct.table_dim.x-2);
    const float y = curand_uniform(&rand_state)*(mstruct.table_dim.y-1);
    const int iy = clamp(__float2int_rd(y), 0, mstruct.table_dim.y-2);
    const int mid = pstruct.material_idx_dev_p[pid];
    const float omega0 = interp2(mstruct.inelastic_dev_p, mstruct.table_pitch, mstruct.table_dim.y+1, ix, 1+iy, mid, x-ix, y-iy);
    float B = -1.0f;
    if(omega0 > 100.0f) {
        const float x = __fdividef(__logf((omega0+10.0f)/mstruct.K_energy_range.x), __logf(mstruct.K_energy_range.y/mstruct.K_energy_range.x))*(mstruct.table_dim.x-1);
        const int ix = clamp(__float2int_rd(x), 0, mstruct.table_dim.x-2);
        const float y = (0.5f+curand_uniform(&rand_state))*(mstruct.table_dim.y-1);
        const int iy = __float2int_rd(y);
        B = cuda_make_ptr(mstruct.ionization_dev_p, mstruct.table_pitch, mstruct.table_dim.y, iy, mid)[ix];
        if(B < 50.0f)
            B = -1.0f;
    }
    if(B < 0) {
        if(mid == 0) {
            /* silicon */
            if(omega0 > 100.0f)
                B = 100.0f;
            else if(omega0 > 8.9f)
                B = 8.9f;
            else if(omega0 > 5.0f)
                B = 5.0f;
            else if(omega0> 1.12f)
                B = 1.12f;
        } else if(mid == 1) {
            /* pmma */
            if(omega0 > 5.0f)
                B = 5.0f;
            else if(omega0 > 3.0f)
                B = 3.0f;
        }
    }

    const float F = mstruct.fermi_dev_p[mid];
    float omega_max = 0.5f*(K+omega0-F); // upper limit of Eq.9 (Ashley), but corrected for the fermi energy.
    float omega_min = omega0;
    float w0 = min(omega0-1.0f, fmaxf(0.0f, B)-F);
    if(K > 2.0f*omega0) {
        omega_min = 0.5f*K*(1.0f-sqrtf(1.0f-2.0f*omega0/K)+omega0/K); // Eq. 10 (Ashley)
        w0 = omega0;
    }
    float omega;
    if((w0 > 0) && (omega_min > w0) && (omega_min < omega_max)) {
        // For nonzero binding energy, sample omega according to eq. 7 in Ashley,
        // using the lower and upper limits as defined above.
        // For inner-shell ionization (Ebind > 50 eV) we substitute the Fermi-energy corrected
        // binding energy for omegaprime (so that the differential cross section becomes inversely
        // proportional to both the total energy transfer and the kinetic energy of the secondary
        // electron).
        const float U = curand_uniform(&rand_state);
        omega = w0/(1.0f-(1.0f-w0/omega_min)*__expf(U*log1pf(-w0/omega_max))*__expf(-U*log1pf(-w0/omega_min)));
    } else {
        // In some cases (typically only occuring for B < 50 eV) we get omega_min > omega_max.
        // This is due to our Fermi energy correction in the definition of omega_max. Physically, this
        // means that momentum cannot be conserved because the primary electron cannot have a final
        // kinetic energy that is lower than the Fermi energy. In this (relatively rare) case we have
        // to ignore momentum conservation and probe omega according to a 1/(omega)^2 distribution
        // with omega0 and omega_max as lower and upper limits, respectively.
        const float U = curand_uniform(&rand_state);
        omega = omega0/(1.0f-U*(1.0f-omega0/omega_max));
    }

    if(B < 0) {
        const float G = mstruct.bandgap_dev_p[mid];
        if(G < 0) {
            // metal: excitation of a fermi sea electron
            /*
            const float z = F/omega;
            if((z-1.0f > fermi_sea::z1) && (z < fermi_sea::z2)) {
                const float x1 = logscale(z-1.0f,fermi_sea::z1,fermi_sea::z2);
                const float x2 = logscale(z,fermi_sea::z1,fermi_sea::z2);
                const float p1 = __expf(interp1(fermi_sea::z2p,fermi_sea::n,x1));
                const float p2 = __expf(interp1(fermi_sea::z2p,fermi_sea::n,x2));
                const float U = curand_uniform(&rand_state);
                const float y = logscale(p1*(1.0f-U)+p2*U,fermi_sea::p1,fermi_sea::p2);
                B = F-omega*__expf(interp1(fermi_sea::p2z,fermi_sea::n,y));
            }
            */
        } else if(omega0 > G) {
            // cross-bandgap excitation
            B = G;
        } else {
            // phonon loss
            pstruct.K_energy_dev_p[pid] = K-omega0;
            rand_state_dev_p[pid] = rand_state;
            return;
        }
    }

    B = fmaxf(0.0f, B);
    const float _K = K-F+2.0f*B;
    const float dK = B+omega;
    const float cos_alpha = __saturatef(sqrtf(__fdividef((1.0f-dK/_K)*(1.0f+0.5f*_K/mc2), 1.0f+0.5f*(_K-dK)/mc2)));
    const float sin_alpha = sqrtf(1.0f-cos_alpha*cos_alpha);
    const float cos_beta = __saturatef(sqrtf(__fdividef((dK/_K)*(1.0f+0.5f*_K/mc2), 1.0f+0.5f*dK/mc2)));
    const float sin_beta = sqrtf(1.0f-cos_beta*cos_beta);
    pstruct.status_dev_p[sid] = cuda_particle_struct::NEW_SECONDARY;
    pstruct.material_idx_dev_p[sid] = mid;
    pstruct.particle_tag_dev_p[sid] = pstruct.particle_tag_dev_p[pid];
    pstruct.pos_x_dev_p[sid] = pstruct.pos_x_dev_p[pid];
    pstruct.pos_y_dev_p[sid] = pstruct.pos_y_dev_p[pid];
    pstruct.pos_z_dev_p[sid] = pstruct.pos_z_dev_p[pid];
    float3 dir = make_float3(pstruct.dir_x_dev_p[pid], pstruct.dir_y_dev_p[pid], pstruct.dir_z_dev_p[pid]);
    dir = dir*rnorm3df(dir.x, dir.y, dir.z);
    float sin_azimuth, cos_azimuth;
    __sincosf(atan2f(dir.y, dir.x), &sin_azimuth, &cos_azimuth);
    const float3 unit_v = make_float3(dir.z*cos_azimuth, dir.z*sin_azimuth, -sqrtf(__saturatef(1.0f-dir.z*dir.z)));
    const float3 unit_u = cross_product(unit_v, dir);
    float sin_phi, cos_phi;
    sincospif(2.0f*curand_uniform(&rand_state), &sin_phi, &cos_phi);
    pstruct.K_energy_dev_p[pid] = K-omega;
    pstruct.dir_x_dev_p[pid] = dir.x*cos_alpha+(unit_u.x*cos_phi+unit_v.x*sin_phi)*sin_alpha;
    pstruct.dir_y_dev_p[pid] = dir.y*cos_alpha+(unit_u.y*cos_phi+unit_v.y*sin_phi)*sin_alpha;
    pstruct.dir_z_dev_p[pid] = dir.z*cos_alpha+(unit_u.z*cos_phi+unit_v.z*sin_phi)*sin_alpha;
    pstruct.K_energy_dev_p[sid] = F+omega-B;
    pstruct.dir_x_dev_p[sid] = dir.x*cos_beta-(unit_u.x*cos_phi+unit_v.x*sin_phi)*sin_beta;
    pstruct.dir_y_dev_p[sid] = dir.y*cos_beta-(unit_u.y*cos_phi+unit_v.y*sin_phi)*sin_beta;
    pstruct.dir_z_dev_p[sid] = dir.z*cos_beta-(unit_u.z*cos_phi+unit_v.z*sin_phi)*sin_beta;
    rand_state_dev_p[pid] = rand_state;
}
