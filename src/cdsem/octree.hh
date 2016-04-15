/**
 * @file src/cdsem/octree.hh
 * @author Thomas Verduin <T.Verduin@tudelft.nl>
 */

#ifndef eSCATTER__CDSEM__OCTREE__HEADER_INCLUDED
#define eSCATTER__CDSEM__OCTREE__HEADER_INCLUDED

#include <cinttypes>
#include <functional>
#include <vector>
#include "point3.hh"
#include "triangle.hh"

class octree {
public:
    using triangle_p_vector = std::vector<const::triangle*>;

    octree(const point3& AABB_min, const point3& AABB_max);
    octree(const octree&);
    octree& operator=(const octree&) = delete;
    ~octree();

    const point3& center() const;
    const point3& halfsize() const;
    bool leaf() const;
    bool empty() const;
    int count() const;

    int level() const;
    int depth() const;
    int occupancy() const;
    uint64_t location() const;

    const triangle* insert(const triangle&);
    std::pair<const triangle*,double> intersect(const point3&, const point3&) const;

    octree* root();
    const octree* root() const;
    const octree* parent() const;
    const octree* traverse(int octant) const;

    triangle_p_vector::const_iterator cbegin() const;
    triangle_p_vector::const_iterator cend() const;

    int octant(const point3&) const;
    bool overlap(const triangle&) const;
    bool overlap(const point3&, const point3&) const;

public:
    static const int _max_count = 16;
    static const int _max_depth = 21;
    point3 _AABB_center;
    point3 _AABB_halfsize;
    octree* _parent_p;
    octree* _child_p[8];
    triangle_p_vector _triangles;
};

#endif