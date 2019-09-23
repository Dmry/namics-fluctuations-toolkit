#include "lattice_accessor.h"

Lattice_accessor::Lattice_accessor() :
MX{0}, MY{0}, MZ{0}, dimensionality{static_cast<Dimensionality>(1)} 
{}

void Lattice_accessor::set_jumps() noexcept {
    switch (dimensionality) {
    case 1:
        jump_x = 1;
        jump_y = 0;
        jump_z = 0;
        break;
    case 2:
        jump_x = (MY+BOUNDARIES);
        jump_y = 1;
        jump_z = 0;
        break;
    case 3:
        jump_x = (MY+BOUNDARIES)*(MZ+BOUNDARIES);
        jump_y = (MZ+BOUNDARIES);
        jump_z = 1;
        break;
    }
    system_size = (MZ+BOUNDARIES)*(MY+BOUNDARIES)*(MX+BOUNDARIES);
}

const Lattice_accessor::Coordinate Lattice_accessor::coordinate(size_t index) {
    
    size_t mod = 0;
    Coordinate coordinate;

    coordinate[Dimension::X] = index / jump_x;
    mod = index % jump_x;

    if (dimensionality > 1) {
        coordinate[Dimension::Y] = mod / jump_y;
        mod = mod % jump_y;
    }

    if (dimensionality > 2) {
        coordinate[Dimension::Z] = mod / jump_z;
    }
    return coordinate;
}

void Lattice_accessor::skip_bounds(std::function<void(size_t, size_t, size_t)> function) noexcept {
size_t x{0};
size_t y{0};

    for ( size_t z = SYSTEM_EDGE_OFFSET ; z < MZ+SYSTEM_EDGE_OFFSET ; ++z ) {
        y = SYSTEM_EDGE_OFFSET;
        do {
            x = SYSTEM_EDGE_OFFSET;
            do {
                function(x, y, z);
                ++x;
            } while (x < MX+SYSTEM_EDGE_OFFSET );
            ++y;
        } while (y < MY+SYSTEM_EDGE_OFFSET );
    }
}

void Lattice_accessor::system_plus_bounds(std::function<void(size_t, size_t, size_t)> function) noexcept {
size_t x{0};
size_t y{0};

    for ( size_t z = 0 ; z < MZ+BOUNDARIES ; ++z ) {
        y = 0;
        do {
            x = 0;
            do {
                function(x, y, z);
                ++x;
            } while (x < MX+BOUNDARIES );
            ++y;
        } while (y < MY+BOUNDARIES );
    }
}

size_t Lattice_accessor::index (const size_t x, const size_t y, const size_t z) noexcept {
    return (x*jump_x + y*jump_y + z*jump_z);
}

size_t Lattice_accessor::index (const size_t x, const size_t y, const size_t z) const noexcept {
    return (x*jump_x + y*jump_y + z*jump_z);
}


void Lattice_accessor::x0_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept {
size_t x = 0, y = 0, z = 0;
do {
    z = 0;
    do {
        function(x,y,z);
        ++z;
        } while (z < MZ+BOUNDARIES);
    ++y;
    } while (y < MY+BOUNDARIES);   
}

void Lattice_accessor::xm_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept {
    size_t x = MX + SYSTEM_EDGE_OFFSET, y = 0, z = 0;
    do {
    z = 0;
    do {
        function(x,y,z);
        ++z;
    } while (z < MZ+BOUNDARIES);
    ++y;
    } while (y < MY+BOUNDARIES);
}

void Lattice_accessor::y0_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept {
    size_t x = 0, y = 0, z = 0;
    do {
        z = 0;
        do {
            function(x,y,z);
            ++z;
        } while (z < MZ+BOUNDARIES);
    ++x;
    } while (x < MX+BOUNDARIES);
}

void Lattice_accessor::ym_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept {
    size_t x = 0, y = MY + SYSTEM_EDGE_OFFSET, z = 0;
    do {
        z = 0;
        do {
            function(x,y,z);
            ++z;
        } while (z < MZ+BOUNDARIES);
        ++x;
    } while (x < MX+BOUNDARIES);
}

void Lattice_accessor::z0_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept {
    size_t x = 0, y = 0, z = 0;
    do {
        y = 0;
        do {
            function(x,y,z);
            ++y;
        } while (y < MY+BOUNDARIES);
        ++x;
    } while (x < MX+BOUNDARIES);
}

void Lattice_accessor::zm_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept {
    size_t x = 0, y = 0, z = MZ + SYSTEM_EDGE_OFFSET;
    do {
        y = 0;
        do {
            function(x,y,z);
            ++y;
        } while (y < MY+BOUNDARIES);
        ++x;
    } while (x < MX+BOUNDARIES);
}