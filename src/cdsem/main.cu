/**
 * @file src/cdsem/main.cu
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 * @author Sebastiaan Lokhorst <S.R.Lokhorst@tudelft.nl>
 */

#include <algorithm>
#include <chrono>
#include <cfloat>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <curand_kernel.h>
#include <cub/cub.cuh>
#include <common/archive.hh>
#include <common/cuda_mem_scope.cuh>
#include <common/cuda_safe_call.cuh>
#include "cuda_kernels.cuh"
#include "material.hh"
#include "trimesh.hh"

const int cuda_warp_size = 32;
const int cuda_block_size = cuda_warp_size*4;

class sc {
public:
	sc() : _cursor_vec(std::vector<char>{ '|', '/', '-', '\\' }) {}
    void reset() {
       _cursor_index = 0;
    }
    void print() {
        std::clog << "    \r [" << _cursor_vec[(_cursor_index++)%_cursor_vec.size()] << "]";
    }
    void finish() const {
        std::clog << "    \r [*]";
    }
private:
    std::vector<char> _cursor_vec;
    size_t _cursor_index = 0;
} spin_cursor;


int main(const int argc, char* argv[]) {
    if(argc < 3)
        return EXIT_FAILURE;

    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
    });

    const int grid_cell_count = 1.6e6;
    const int electron_capacity = 2e5;
    const int cuda_block_count = 1+electron_capacity/cuda_block_size;
    const int prescan_size = 256;

    trimesh triangle_mesh;
    std::ifstream tri_ifs(argv[1]);
    std::clog << " [*] loading triangles";
    std::clog << " file='" << argv[1] << "'";
    std::clog << std::endl;
    while(!tri_ifs.eof()) {
        trimesh::face triangle;
        tri_ifs >> triangle.in >> triangle.out;
        if(tri_ifs.eof())
            break;
        tri_ifs >> triangle.A.x >> triangle.A.y >> triangle.A.z;
        tri_ifs >> triangle.B.x >> triangle.B.y >> triangle.B.z;
        tri_ifs >> triangle.C.x >> triangle.C.y >> triangle.C.z;
        triangle_mesh.push(triangle);
    }
    tri_ifs.close();

    std::clog << " [*] initializing geometry";
    cuda_geometry_struct gstruct(cuda_geometry_struct::create(triangle_mesh, grid_cell_count));
    std::clog << " cell=" << gstruct.cell.x << "x" << gstruct.cell.y << "x" << gstruct.cell.z;
    std::clog << " dim=" << gstruct.dim.x << "x" << gstruct.dim.y << "x" << gstruct.dim.z << "x" << gstruct.dim.w;
    std::clog << std::endl;

    std::vector<material> mat_vec;
    for(const std::string& mat_file : {"../data/silicon.mat", "../data/pmma.mat"}) {
        std::clog << " [*] loading material";
        std::clog << " index=" << mat_vec.size();
        std::clog << " file='" << mat_file << "'";
        std::clog << std::endl;
        std::ifstream ifs(mat_file, std::ifstream::in | std::ifstream::binary);
        if(!ifs.is_open())
            throw std::ios_base::failure("failed to open '"+mat_file+"' for reading");
        archive::istream ia(ifs);
        mat_vec.push_back(material());
        ia >> mat_vec.back();
        ifs.close();
    }

    cuda_material_struct mstruct(cuda_material_struct::create(mat_vec.size()));
    for(int i = 0; i < mstruct.n; i++)
        mstruct.assign(i, mat_vec[i]);

    cuda_particle_struct pstruct(cuda_particle_struct::create(electron_capacity));
    struct primary {
        float3 pos;
        float3 dir;
        float K;
        int tag;
    };
    std::vector<primary> primary_vec;
    std::map<int,std::pair<int,int>> tag_map;
    std::ifstream pri_ifs(argv[2]);
    spin_cursor.reset();
    while(!pri_ifs.eof()) {
        primary pe;
        pri_ifs.read(reinterpret_cast<char*>(&(pe.pos.x)), sizeof(pe.pos.x));
        pri_ifs.read(reinterpret_cast<char*>(&(pe.pos.y)), sizeof(pe.pos.y));
        pri_ifs.read(reinterpret_cast<char*>(&(pe.pos.z)), sizeof(pe.pos.z));
        pri_ifs.read(reinterpret_cast<char*>(&(pe.dir.x)), sizeof(pe.dir.x));
        pri_ifs.read(reinterpret_cast<char*>(&(pe.dir.y)), sizeof(pe.dir.y));
        pri_ifs.read(reinterpret_cast<char*>(&(pe.dir.z)), sizeof(pe.dir.z));
        pri_ifs.read(reinterpret_cast<char*>(&pe.K), sizeof(pe.K));
        int2 pixel;
        pri_ifs.read(reinterpret_cast<char*>(&(pixel.x)), sizeof(pixel.x));
        pri_ifs.read(reinterpret_cast<char*>(&(pixel.y)), sizeof(pixel.y));
        pe.tag = tag_map.size();
        if((pe.tag%65536 == 0) || (pri_ifs.eof())) {
            spin_cursor.print();
            std::clog << " loading primary electrons";
            std::clog << " file='" << argv[2] << "'";
            std::clog << " count=" << primary_vec.size();
            std::clog << std::flush;
        }
        if(pri_ifs.eof())
            break;
        if(pe.pos.x < -32)
            pe.pos.x = -64-pe.pos.x;
        if(pe.pos.x > 32)
            pe.pos.x = 64-pe.pos.x;
        if(pe.pos.y < -750)
            pe.pos.y = -1500-pe.pos.y;
        if(pe.pos.y > 750)
            pe.pos.y = 1500-pe.pos.y;
        primary_vec.push_back(pe);
        tag_map[pe.tag] = std::make_pair(pixel.x, pixel.y);
    }
    spin_cursor.finish();
    std::clog << std::endl;
    pri_ifs.close();
    const int np = primary_vec.size();
    std::vector<float3> pos_vec;
    std::vector<float3> dir_vec;
    std::vector<float> K_vec;
    std::vector<int> tag_vec;
    const unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::shuffle(primary_vec.begin(), primary_vec.end(), std::default_random_engine(seed));
    int i;
    for(i = 0; i < prescan_size; i++) {
        pos_vec.push_back(primary_vec[i].pos);
        dir_vec.push_back(primary_vec[i].dir);
        K_vec.push_back(primary_vec[i].K);
        tag_vec.push_back(primary_vec[i].tag);
    }
    std::multimap<int,int> primary_map;
    spin_cursor.reset();
    for(; i < np; i++) {
        int4 gid;
        gid.x = std::floor((primary_vec[i].pos.x-gstruct.org.x)/gstruct.cell.x);
        gid.y = std::floor((primary_vec[i].pos.y-gstruct.org.y)/gstruct.cell.y);
        gid.z = std::floor((primary_vec[i].pos.z-gstruct.org.z)/gstruct.cell.z);
        gid.w = gid.x+gstruct.dim.x*gid.y+gstruct.dim.x*gstruct.dim.y*gid.z;
        primary_map.insert(std::make_pair(gid.w, i));
        if((i%65536 == 0) || (i == np-1)) {
            spin_cursor.print();
            std::clog << " sorting primary electrons";
            std::clog << " pct=" << 100*i/(np-1) << "%";
            std::clog << std::flush;
        }
    }
    spin_cursor.finish();
    std::clog << std::endl;
    for(auto cit = primary_map.cbegin(); cit != primary_map.cend(); cit++) {
        const int i = cit->second;
        pos_vec.push_back(primary_vec[i].pos);
        dir_vec.push_back(primary_vec[i].dir);
        K_vec.push_back(primary_vec[i].K);
        tag_vec.push_back(primary_vec[i].tag);
    }

    int* radix_index_dev_p = nullptr;
    void* radix_dump_dev_p = nullptr;
    void* radix_temp_dev_p = nullptr;
    size_t radix_temp_size = 0;
    curandState* rand_state_dev_p = nullptr;
    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaMalloc(&radix_index_dev_p, electron_capacity*sizeof(int));
        cuda_mem_scope<int>(radix_index_dev_p, electron_capacity, [&](int* index_p) {
            for(int i = 0; i < electron_capacity; i++)
                index_p[i] = i;
        });
        cudaMalloc(&radix_dump_dev_p, electron_capacity*sizeof(int));
        cub::DeviceRadixSort::SortPairs<int,int>(
            nullptr, radix_temp_size,
            pstruct.status_dev_p, static_cast<int*>(radix_dump_dev_p),
            radix_index_dev_p, pstruct.pid_dev_p,
            electron_capacity
        );
        radix_temp_size++;
        cudaMalloc(&radix_temp_dev_p, radix_temp_size);
        cudaMalloc(&rand_state_dev_p, electron_capacity*sizeof(curandState));
        __init_rand_state<<<cuda_block_count,cuda_block_size>>>(
            rand_state_dev_p, 0, electron_capacity
        );
    });

    auto __execute_iteration = [&] {
        cuda_safe_call(__FILE__, __LINE__, [&]() {
            __init_trajectory<<<cuda_block_count,cuda_block_size>>>(pstruct, gstruct, mstruct);
            __probe_isec_event<<<cuda_block_count,cuda_block_size>>>(pstruct, gstruct);
            __probe_scatter_event<<<cuda_block_count,cuda_block_size>>>(pstruct, mstruct, rand_state_dev_p);
            __update_trajectory<<<cuda_block_count,cuda_block_size>>>(pstruct, gstruct);
            cub::DeviceRadixSort::SortPairs<int,int>(
                radix_temp_dev_p, radix_temp_size,
                pstruct.status_dev_p, static_cast<int*>(radix_dump_dev_p),
                radix_index_dev_p, pstruct.pid_dev_p,
                electron_capacity, 0, 4
            );
            __apply_isec_event<<<cuda_block_count,cuda_block_size>>>(pstruct, gstruct, mstruct, rand_state_dev_p);
            __apply_elastic_event<<<cuda_block_count,cuda_block_size>>>(pstruct, mstruct, rand_state_dev_p);
            __apply_inelastic_event<<<cuda_block_count,cuda_block_size>>>(pstruct, mstruct, rand_state_dev_p);
        });
    };

    auto __flush_detected_particles = [&](std::ostream& os) {
        pstruct.for_each(cuda_particle_struct::DETECTED, [&](float3 pos, float3 dir, float K, int tag) {
            const int2 pixel = make_int2(tag_map[tag].first, tag_map[tag].second);
            os.write(reinterpret_cast<const char*>(&(pos.x)), sizeof(pos.x));
            os.write(reinterpret_cast<const char*>(&(pos.y)), sizeof(pos.y));
            os.write(reinterpret_cast<const char*>(&(pos.z)), sizeof(pos.z));
            os.write(reinterpret_cast<const char*>(&(dir.x)), sizeof(dir.x));
            os.write(reinterpret_cast<const char*>(&(dir.y)), sizeof(dir.y));
            os.write(reinterpret_cast<const char*>(&(dir.z)), sizeof(dir.z));
            os.write(reinterpret_cast<const char*>(&K), sizeof(K));
            os.write(reinterpret_cast<const char*>(&(pixel.x)), sizeof(pixel.x));
            os.write(reinterpret_cast<const char*>(&(pixel.y)), sizeof(pixel.y));
        });
        pstruct.flush();
    };

    int batch_index = pstruct.push(pos_vec.data(), dir_vec.data(), K_vec.data(), tag_vec.data(), prescan_size);
    std::vector<std::pair<int,int>> prescan_stats_vec;
    prescan_stats_vec.push_back(std::make_pair(batch_index, 0));
    spin_cursor.reset();
    while(prescan_stats_vec.back().first > 0) {
        __execute_iteration();
        cuda_safe_call(__FILE__, __LINE__, [&]() {
            cuda_mem_scope<int>(pstruct.status_dev_p, electron_capacity, [&](int* status_p) {
                int running_count = 0;
                int detected_count = 0;
                for(int i = 0; i < electron_capacity; i++)
                    switch(status_p[i]) {
                        case cuda_particle_struct::DETECTED:
                            detected_count++;
                        break;
                        case cuda_particle_struct::TERMINATED:
                        break;
                        default:
                            running_count++;
                        break;
                    }
                prescan_stats_vec.push_back(std::make_pair(running_count, detected_count));
            }); // status_dev_p
        });
        spin_cursor.print();
        std::clog << " executing prescan";
        std::clog << " running_count=" << prescan_stats_vec.back().first;
        std::clog << " detected_count=" << prescan_stats_vec.back().second;
        std::clog << std::flush;
    }
    spin_cursor.finish();
    std::clog << std::endl;
    int frame_size = 0;
    for(size_t i = 0; i < prescan_stats_vec.size(); i++)
        if(prescan_stats_vec[i].first > prescan_stats_vec[frame_size].first)
            frame_size = i;
    float accumulator = 0;
    for(size_t i = 2*frame_size; i < prescan_stats_vec.size(); i += frame_size)
        accumulator += 1.0*prescan_stats_vec[i].first/prescan_size;
    accumulator += 2.0*prescan_stats_vec[frame_size].first/prescan_size;
    accumulator += 2.0*prescan_stats_vec[frame_size].second/prescan_size;
    int batch_size = 0.95*electron_capacity/accumulator;

    const std::chrono::high_resolution_clock::time_point chrono_start = std::chrono::high_resolution_clock::now();
    int running_count;
    spin_cursor.reset();
    do {
        batch_size = std::min(np-batch_index, batch_size);
        if(batch_size > 0)
            batch_index += pstruct.push(pos_vec.data()+batch_index, dir_vec.data()+batch_index, K_vec.data()+batch_index, tag_vec.data()+batch_index, batch_size);
        for(int i = 0; i < frame_size; i++)
            __execute_iteration();
        __flush_detected_particles(std::cout);
        running_count = 0;
        cuda_safe_call(__FILE__, __LINE__, [&]() {
            cuda_mem_scope<int>(pstruct.status_dev_p, electron_capacity, [&](int* status_p) {
                for(int i = 0; i < electron_capacity; i++)
                    switch(status_p[i]) {
                        case cuda_particle_struct::TERMINATED:
                        break;
                        default:
                            running_count++;
                        break;
                    }
            }); // status_dev_p
        });
        spin_cursor.print();
        std::clog << " executing exposure";
        std::clog << " pct=" << 100*batch_index/(np-1) << "%";
        std::clog << " frame_size=" << frame_size;
        std::clog << " batch_size=" << batch_size;
        std::clog << " running_count=" << running_count;
        std::clog << std::flush;
    } while(running_count > 0);
    spin_cursor.finish();
    std::clog << std::endl;
    const std::chrono::high_resolution_clock::time_point chrono_stop = std::chrono::high_resolution_clock::now();
    std::clog << " [*] time_lapse=" << std::chrono::duration_cast<std::chrono::seconds>(chrono_stop-chrono_start).count() << " seconds";
    std::clog << std::endl;

    cuda_safe_call(__FILE__, __LINE__, [&]() {
        cudaFree(radix_temp_dev_p);
        cudaFree(radix_dump_dev_p);
        cudaFree(radix_index_dev_p);
        cudaFree(rand_state_dev_p);
    });
    cuda_geometry_struct::release(gstruct);
    cuda_material_struct::release(mstruct);
    cuda_particle_struct::release(pstruct);

    return EXIT_SUCCESS;
}