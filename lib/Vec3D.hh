// Vec3D.hh
// T. M. Kelley
// Dec 19, 2012
// (c) Copyright 2012 LANSLLC, all rights reserved


#ifndef VEC3D_HH
#define VEC3D_HH


#include <array>
#include <stdint.h>
#include <iostream> // for operator >>, <<, istream, ostream.

namespace nut
{
    template <typename fp_t, size_t dim>
    struct Vec_T
    {
        typedef fp_t const & gcr;
        typedef fp_t value_type;
        typedef fp_t * iterator;
        typedef fp_t const * const_iterator;

        typedef Vec_T<fp_t, dim> vec_t;

        iterator begin(){return &v[0];}
        iterator end(){return &v[NCS];}
        const_iterator begin() const {return &v[0];}
        const_iterator end() const {return &v[NCS];}

        static const uint32_t NCS=dim;

        std::array<fp_t,dim> v;

        size_t size() const {return size_t(3);}

        // Vec_T(gcr v1,gcr v2, gcr v3): v{v1,v2,v3}{}

        Vec_T(){std::fill(v.begin(),v.end(),fp_t(0.0));}

        Vec_T(fp_t const i) {  std::fill(v.begin(),v.end(),i);}

        Vec_T div_by(fp_t const div) const
        {
            vec_t vout;
            for(uint32_t d = 0; d < dim; ++d)
            {
                vout.v[d] = v[d]/div;
            }
            return std::move(vout);
        }

        // for STL
        vec_t & operator=(vec_t const & rhs)
        {
            if(this == &rhs) return *this;
            std::copy(rhs.begin(),rhs.end(),v.begin());
            return *this;
        }

        bool operator==(vec_t const & rhs) const
        {
            return std::equal(rhs.begin(),rhs.end(),v.begin());
        }

        bool operator!=(vec_t const & rhs) const
        {
            return !(*this == rhs);
        }

        Vec_T & operator+=(vec_t const & rhs)
        {
            for(uint32_t d = 0; d < dim; ++d)
            {
                v[d]+=rhs.v[d];
            }
            return *this;
        }

        Vec_T operator-(vec_t const & rhs)
        {
            for(uint32_t d = 0; d < dim; ++d)
            {
                v[d]-=rhs.v[d];
            }
            return *this;
        }


        friend std::ostream &
        operator<<(std::ostream & s,Vec_T const &v)
        {
            s << "{";
            for(uint32_t d = 0; d < (dim-1); ++d)
            {
                s << v.v[d] << ",";
            }
            s << v.v[dim-1] << "}";
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
                for(uint32_t d = 0; d < dim; ++d)
                {
                    s >> v.v[d] >> c;
                }
            }
            else
            {
                s.putback(c);
                for(uint32_t d = 0; d < dim; ++d)
                {
                    s >> v.v[d];
                }
            }
            return s;
        }
    }; // Vec_T


    /** scale a vector by a scalar, generating a new vector */
    template <typename fp_t,size_t dim>
    inline
    Vec_T<fp_t,dim> operator*(Vec_T<fp_t,dim> const & v, fp_t const f)
    {
        Vec_T<fp_t,dim> vout;
        for(uint32_t d = 0; d < dim; ++d)
        {
            vout.v[d] = v.v[d] * f;
        }
        return std::move(vout);
    }
    /** scale a vector by a scalar, generating a new vector */
    template <typename fp_t,size_t dim>
    inline
    Vec_T<fp_t,dim> operator*(fp_t const f, Vec_T<fp_t,dim> const & v)
    {
        Vec_T<fp_t,dim> vout;
        for(uint32_t d = 0; d < dim; ++d)
        {
            vout.v[d] = v.v[d] * f;
        }
        return std::move(vout);
    }
    /** generate a new vector by unary negating a vector */
    template <typename fp_t,size_t dim>
    inline
    Vec_T<fp_t,dim> operator-(Vec_T<fp_t,dim> const & v)
    {
        Vec_T<fp_t,dim> vout;
        for(uint32_t d = 0; d < dim; ++d)
        {
            vout.v[d] = -v.v[d];
        }
        return std::move(vout);
    }
    /** generate a new vector by adding two vectors */
    template <typename fp_t,size_t dim>
    inline
    Vec_T<fp_t,dim> operator+(Vec_T<fp_t,dim> const & v1,Vec_T<fp_t,dim> const & v2)
    {
        Vec_T<fp_t,dim> vout;
        for(uint32_t d = 0; d < dim; ++d)
        {
            vout.v[d] = v1.v[d] + v2.v[d];
        }
        return std::move(vout);
    }

    /** generate a new vector by subtracting v2 from v1 */
    template <typename fp_t,size_t dim>
    inline
    Vec_T<fp_t,dim> operator-(Vec_T<fp_t,dim> const & v1,Vec_T<fp_t,dim> const & v2)
    {
        Vec_T<fp_t,dim> vout;
        for(uint32_t d = 0; d < dim; ++d)
        {
            vout.v[d] = v1.v[d] - v2.v[d];
        }
        return std::move(vout);
    }

    template <typename fp_t,size_t dim>
    fp_t inline
    dot(Vec_T<fp_t,dim> const & v1, Vec_T<fp_t,dim> const & v2)
    {
        fp_t init(0.0);
        return std::inner_product(v1.begin(),v1.end(),v2.begin(),init);
    }



} // nut::

#endif // include guard


// End of file
