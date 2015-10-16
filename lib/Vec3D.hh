// Vec3D.hh
// T. M. Kelley
// Dec 19, 2012
// (c) Copyright 2012 LANSLLC, all rights reserved


#ifndef VEC3D_HH
#define VEC3D_HH

#include <stdint.h>
#include <iostream> // for operator >>, <<, istream, ostream.

namespace nut
{
    template <typename fp_t>
    struct Vec_T
    {
        typedef fp_t const & gcr;
        typedef fp_t value_type;
        typedef fp_t * iterator;
        typedef fp_t const * const_iterator;
        iterator begin(){return &v[0];}
        iterator end(){return &v[NCS];}
        const_iterator begin() const {return &v[0];}
        const_iterator end() const {return &v[NCS];}

        static const uint32_t NCS=3;

        fp_t v[3];

        size_t size() const {return size_t(3);}

        Vec_T(gcr v1,gcr v2, gcr v3): v{v1,v2,v3}{}

        Vec_T(): v{0,0,0}{}

        Vec_T(fp_t i) : v{i,i,i}{}

        Vec_T div_by(fp_t const d) const
        {
            return std::move(Vec_T{v[0]/d,v[1]/d,v[2]/d});
        }

        // for STL
        Vec_T & operator=(Vec_T const & rhs)
        {
            if(this == &rhs) return *this;
            v[0]=rhs.v[0];
            v[1]=rhs.v[1];
            v[2]=rhs.v[2];
            return *this;
        }

        bool operator==(Vec_T const & rhs) const
        {
            return v[0] == rhs.v[0] && v[1] == rhs.v[1] && v[2] == rhs.v[2];
        }

        bool operator!=(Vec_T const & rhs) const
        {
            return !(*this == rhs);
        }

        Vec_T & operator+=(Vec_T const & rhs)
        {
            v[0]+=rhs.v[0];
            v[1]+=rhs.v[1];
            v[2]+=rhs.v[2];
            return *this;
        }

        Vec_T operator-(Vec_T const & rhs)
        {
            return std::move(Vec_T{v[0]-rhs.v[0],v[1]-rhs.v[1],v[2]-rhs.v[2]});
        }


        friend std::ostream &
        operator<<(std::ostream & s,Vec_T const &v)
        {
            s << "{" << v.v[0] << "," << v.v[1] << "," << v.v[2] << "}";
            return s;
        }

        friend std::istream &
        operator>>(std::istream & s,Vec_T & v)
        {
            // read two possible formats
            // <w>vx<w>vy<w>vz
            //       or
            // {vx,vy,vz}
            char c;
            s >> c;
            if(c == '{')
            {
                s >> v.v[0] >> c;
                s >> v.v[1] >> c;
                s >> v.v[2] >> c;
            }
            else
            {
                s.putback(c);
                s >> v.v[0] >> v.v[1] >> v.v[2];
            }
            return s;
        }
    }; // Vec_T

    /** scale a vector by a scalar, generating a new vector */
    template <typename fp_t>
    Vec_T<fp_t> operator*(Vec_T<fp_t> const & v, fp_t const f)
    {
        return std::move(Vec_T<fp_t>{v.v[0]*f,v.v[1]*f,v.v[2]*f});
    }
    /** scale a vector by a scalar, generating a new vector */
    template <typename fp_t>
    Vec_T<fp_t> operator*(fp_t const f, Vec_T<fp_t> const & v)
    {
        return std::move(Vec_T<fp_t>{v.v[0]*f,v.v[1]*f,v.v[2]*f});
    }
    /** generate a new vector by unary negating a vector */
    template <typename fp_t>
    Vec_T<fp_t> operator-(Vec_T<fp_t> const & v)
    {
        return std::move(Vec_T<fp_t>{-v.v[0],-v.v[1],-v.v[2]});
    }
    /** generate a new vector by adding two vectors */
    template <typename fp_t>
    Vec_T<fp_t> operator+(Vec_T<fp_t> const & v,Vec_T<fp_t> const & rhs)
    {
        return std::move(
            Vec_T<fp_t>{ v.v[0]+rhs.v[0], v.v[1]+rhs.v[1], v.v[2]+rhs.v[2] }
            );
    }
} // nut::

#endif // include guard


// End of file
