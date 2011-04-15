/*
 *  This file is part of Healpix_cxx.
 *
 *  Healpix_cxx is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  Healpix_cxx is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Healpix_cxx; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix, see http://healpix.jpl.nasa.gov
 */

/*
 *  Healpix_cxx is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*
 *  Copyright (C) 2003-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#include "healpix_base.h"
#include "geom_utils.h"
#include "lsconstants.h"

using namespace std;

int Healpix_Base::nside2order (int nside)
  {
  planck_assert (nside>0, "invalid value for Nside");
  return ((nside)&(nside-1)) ? -1 : ilog2(nside);
  }
int Healpix_Base::npix2nside (int npix)
  {
  int res=isqrt(npix/12);
  planck_assert (npix==res*res*12, "invalid value for npix");
  return res;
  }

int Healpix_Base::ring_above (double z) const
  {
  double az=abs(z);
  if (az<=twothird) // equatorial region
    return int(nside_*(2-1.5*z));
  int iring = int(nside_*sqrt(3*(1-az)));
  return (z>0) ? iring : 4*nside_-iring-1;
  }

void Healpix_Base::in_ring(int iz, double phi0, double dphi,
  rangeset<int> &pixset) const
  {
  int nr, ipix1;
  bool shifted;

  get_ring_info_small(iz,ipix1,nr,shifted);
  double shift = shifted ? 0.5: 0.;

  int ipix2 = ipix1 + nr - 1;       // highest pixel number in the ring

  if (dphi > (pi-1e-12))
    pixset.append(ipix1,ipix2+1);
  else
    {
    int ip_lo = ifloor<int>(nr*inv_twopi*(phi0-dphi) - shift)+1;
    int ip_hi = ifloor<int>(nr*inv_twopi*(phi0+dphi) - shift);
    if (ip_lo<0)
      {
      pixset.append(ipix1,ipix1+ip_hi+1);
      pixset.append(ipix1+ip_lo+nr,ipix2+1);
      }
    else if (ip_hi>=nr)
      {
      pixset.append(ipix1,ipix1+ip_hi-nr+1);
      pixset.append(ipix1+ip_lo,ipix2+1);
      }
    else
      pixset.append(ipix1+ip_lo,ipix1+ip_hi+1);
    }
  }

namespace {

/* Short note on the "zone":
   zone = 0: pixel lies completely outside the queried shape
          1: pixel may overlap with the shape, pixel center is outside
          2: pixel center is inside the shape, but maybe not the complete pixel
          3: pixel lies completely inside the shape */

inline void check_pixel (int o, int order_, int omax, int zone,
  rangeset<int> &pixset, int pix, vector<pair<int,int> > &stk, bool inclusive,
  int &stacktop)
  {
  if (zone==0) return;

  if (o<order_)
    {
    if (zone>=3)
      {
      int sdist=2*(order_-o); // the "bit-shift distance" between map orders
      pixset.append(pix<<sdist,(pix+1)<<sdist); // output all subpixels
      }
    else // (zone>=1)
      for (int i=0; i<4; ++i)
        stk.push_back(make_pair(4*pix+3-i,o+1)); // add children
    }
  else if (o>order_) // this implies that inclusive==true
    {
    if (zone>=2) // pixel center in shape
      {
      pixset.append(pix>>(2*(o-order_))); // output the parent pixel at order_
      stk.resize(stacktop); // unwind the stack
      }
    else // (zone>=1): pixel center in safety range
      {
      if (o<omax) // check sublevels
        for (int i=0; i<4; ++i) // add children in reverse order
          stk.push_back(make_pair(4*pix+3-i,o+1));
      else // at resolution limit
        {
        pixset.append(pix>>(2*(o-order_))); // output the parent pixel at order_
        stk.resize(stacktop); // unwind the stack
        }
      }
    }
  else // o==order_
    {
    if (zone>=2)
      pixset.append(pix);
    else if (inclusive) // and (zone>=1)
      {
      if (order_<omax) // check sublevels
        {
        stacktop=stk.size(); // remember current stack position
        for (int i=0; i<4; ++i) // add children in reverse order
          stk.push_back(make_pair(4*pix+3-i,o+1));
        }
      else // at resolution limit
        pixset.append(pix); // output the pixel
      }
    }
  }

} // unnamed namespace

void Healpix_Base::query_disc (pointing ptg, double radius, bool inclusive,
  rangeset<int> &pixset) const
  {
  pixset.clear();
  ptg.normalize();

  if (scheme_==RING)
    {
    if (inclusive) radius+=max_pixrad();
    if (radius>=pi)
      { pixset.append(0,npix_); return; }

    double cosang = cos(radius);

    double z0 = cos(ptg.theta);
    double xa = 1./sqrt((1-z0)*(1+z0));

    double rlat1  = ptg.theta - radius;
    double zmax = cos(rlat1);
    int irmin = ring_above (zmax)+1;

    if ((rlat1<=0) && (irmin>1)) // north pole in the disk
      {
      int sp,rp; bool dummy;
      get_ring_info_small(irmin-1,sp,rp,dummy);
      pixset.append(0,sp+rp);
      }

    double rlat2  = ptg.theta + radius;
    double zmin = cos(rlat2);
    int irmax = ring_above (zmin);

    for (int iz=irmin; iz<=irmax; ++iz) // rings partially in the disk
      {
      double z=ring2z(iz);

      double x = (cosang-z*z0)*xa;
      double ysq = 1-z*z-x*x;
      planck_assert(ysq>=0, "error in query_disc()");
      double dphi=atan2(sqrt(ysq),x);
      in_ring (iz, ptg.phi, dphi, pixset);
      }

    if ((rlat2>=pi) && (irmax+1<4*nside_)) // south pole in the disk
      {
      int sp,rp; bool dummy;
      get_ring_info_small(irmax+1,sp,rp,dummy);
      pixset.append(sp,npix_);
      }
    }
  else // scheme_==NEST
    {
    if (radius>=pi) // disk covers the whole sphere
      { pixset.append(0,npix_); return; }

    int oplus=inclusive ? 2 : 0;
    int omax=min(int(order_max),order_+oplus); // the order up to which we test

    vec3 vptg(ptg);
    arr<Healpix_Base> base(omax+1);
    arr<double> crpdr(omax+1), crmdr(omax+1);
    double cosrad=cos(radius);
    for (int o=0; o<=omax; ++o) // prepare data at the required orders
      {
      base[o].Set(o,NEST);
      double dr=base[o].max_pixrad(); // safety distance
      crpdr[o] = (radius+dr>pi) ? -1. : cos(radius+dr);
      crmdr[o] = (radius-dr<0.) ?  1. : cos(radius-dr);
      }
    vector<pair<int,int> > stk; // stack for pixel numbers and their orders
    stk.reserve(12+3*omax); // reserve maximum size to avoid reallocation
    for (int i=0; i<12; ++i) // insert the 12 base pixels in reverse order
      stk.push_back(make_pair(11-i,0));

    int stacktop=0; // a place to save a stack position

    while (!stk.empty()) // as long as there are pixels on the stack
      {
      // pop current pixel number and order from the stack
      int pix=stk.back().first,
          o  =stk.back().second;
      stk.pop_back();

      double z,phi;
      base[o].pix2zphi(pix,z,phi);
      // cosine of angular distance between pixel center and disk center
      double cangdist=cosdist_zphi(vptg.z,ptg.phi,z,phi);

      if (cangdist>crpdr[o])
        {
        int zone = (cangdist<cosrad) ? 1 : ((cangdist<=crmdr[o]) ? 2 : 3);

        check_pixel (o, order_, omax, zone, pixset, pix, stk, inclusive,
          stacktop);
        }
      }
    }
  }

void Healpix_Base::query_multidisc (const arr<vec3> &norm,
  const arr<double> &rad, bool inclusive, rangeset<int> &pixset) const
  {
  tsize nv=norm.size();
  planck_assert(nv==rad.size(),"inconsistent input arrays");
  pixset.clear();

  if (scheme_==RING)
    {
    double rplus = inclusive ? max_pixrad() : 0;
    int irmin=1, irmax=4*nside_-1;
    double thmin=0, thmax=pi;
    vector<double> z0,xa,cosrad;
    vector<pointing> ptg;
    z0.reserve(nv); xa.reserve(nv); cosrad.reserve(nv); ptg.reserve(nv);
    for (tsize i=0; i<nv; ++i)
      {
      double r=rad[i]+rplus;
      if (r<pi)
        {
        pointing pnt=pointing(norm[i]);
        cosrad.push_back(cos(r));
        double cth=cos(pnt.theta);
        z0.push_back(cth);
        xa.push_back(1./sqrt((1-cth)*(1+cth)));
        ptg.push_back(pnt);
        double tmp = min(pnt.theta+r,pi);
        if (tmp < thmax)
          { thmax=tmp; irmax=ring_above(cos(thmax)); }
        tmp = max(0.,pnt.theta-r);
        if (tmp > thmin)
          { thmin=tmp; irmin=ring_above(cos(thmin))+1; }
        }
      }

    for (int iz=irmin; iz<=irmax; ++iz)
      {
      double z=ring2z(iz);
      int ipix1,nr;
      bool shifted;
      get_ring_info_small(iz,ipix1,nr,shifted);
      double shift = shifted ? 0.5: 0.;
      int ipix2 = ipix1 + nr - 1;       // highest pixel number in the ring
      rangeset<int> tr;
      tr.append(ipix1,ipix1+nr);
      for (tsize j=0; j<z0.size(); ++j)
        {
        double x = (cosrad[j]-z*z0[j])*xa[j];
        double ysq = 1.-z*z-x*x;
        if (ysq>=0.)
          {
          double dphi=atan2(sqrt(ysq),x);

          if (dphi < (pi-1e-12))
            {
            int ip_lo = ifloor<int>(nr*inv_twopi*(ptg[j].phi-dphi) - shift)+1;
            int ip_hi = ifloor<int>(nr*inv_twopi*(ptg[j].phi+dphi) - shift);
            if (ip_lo<0)
              tr.remove(ipix1+ip_hi+1,ipix1+ip_lo+nr);
            else if (ip_hi>=nr)
              tr.remove(ipix1+ip_hi-nr+1,ipix1+ip_lo);
            else
              {
              tr.remove(ipix1,ipix1+ip_lo);
              tr.remove(ipix1+ip_hi+1,ipix2+1);
              }
            }
          }
        }
      for (tsize j=0; j<tr.size(); ++j)
        pixset.append(tr[j]);
      }
    }
  else // scheme_ == NEST
    {
    int oplus=inclusive ? 2 : 0;
    int omax=min(int(order_max),order_+oplus); // the order up to which we test

    // TODO: ignore all disks with radius>=pi

    arr<Healpix_Base> base(omax+1);
    arr3<double> crlimit(omax+1,nv,3);
    for (int o=0; o<=omax; ++o) // prepare data at the required orders
      {
      base[o].Set(o,NEST);
      double dr=base[o].max_pixrad(); // safety distance
      for (tsize i=0; i<nv; ++i)
        {
        crlimit(o,i,0) = (rad[i]+dr>pi) ? -1. : cos(rad[i]+dr);
        crlimit(o,i,1) = (o==0) ? cos(rad[i]) : crlimit(0,i,1);
        crlimit(o,i,2) = (rad[i]-dr<0.) ?  1. : cos(rad[i]-dr);
        }
      }

    vector<pair<int,int> > stk; // stack for pixel numbers and their orders
    stk.reserve(12+3*omax); // reserve maximum size to avoid reallocation
    for (int i=0; i<12; ++i) // insert the 12 base pixels in reverse order
      stk.push_back(make_pair(11-i,0));

    int stacktop=0; // a place to save a stack position

    while (!stk.empty()) // as long as there are pixels on the stack
      {
      // pop current pixel number and order from the stack
      int pix=stk.back().first,
          o  =stk.back().second;
      stk.pop_back();

      vec3 pv(base[o].pix2vec(pix));

      tsize zone=3;
      for (tsize i=0; i<nv; ++i)
        {
        double crad=dotprod(pv,norm[i]);
        for (tsize iz=0; iz<zone; ++iz)
          if (crad<crlimit(o,i,iz))
            if ((zone=iz)==0) goto bailout;
        }

      check_pixel (o, order_, omax, zone, pixset, pix, stk, inclusive,
        stacktop);
      bailout:;
      }
    }
  }

int Healpix_Base::spread_bits (int v) const
  { return utab[v&0xff] | (utab[(v>>8)&0xff]<<16); }

int Healpix_Base::compress_bits (int v) const
  {
  int raw = (v&0x5555) | ((v&0x55550000)>>15);
  return ctab[raw&0xff] | (ctab[raw>>8]<<4);
  }

void Healpix_Base::nest2xyf (int pix, int &ix, int &iy, int &face_num) const
  {
  face_num = pix>>(2*order_);
  pix &= (npface_-1);
  ix = compress_bits(pix);
  iy = compress_bits(pix>>1);
  }

int Healpix_Base::xyf2nest (int ix, int iy, int face_num) const
  { return (face_num<<(2*order_)) + spread_bits(ix) + (spread_bits(iy)<<1); }

void Healpix_Base::ring2xyf (int pix, int &ix, int &iy, int &face_num) const
  {
  int iring, iphi, kshift, nr;

  int nl2 = 2*nside_;

  if (pix<ncap_) // North Polar cap
    {
    iring = (1+isqrt(1+2*pix))>>1; //counted from North pole
    iphi  = (pix+1) - 2*iring*(iring-1);
    kshift = 0;
    nr = iring;
    face_num=(iphi-1)/nr;
    }
  else if (pix<(npix_-ncap_)) // Equatorial region
    {
    int ip = pix - ncap_;
    int tmp = (order_>=0) ? ip>>(order_+2) : ip/(4*nside_);
    iring = tmp+nside_;
    iphi = ip-tmp*4*nside_ + 1;
    kshift = (iring+nside_)&1;
    nr = nside_;
    unsigned int ire = iring-nside_+1,
                 irm = nl2+2-ire;
    int ifm = iphi - ire/2 + nside_ -1,
        ifp = iphi - irm/2 + nside_ -1;
    if (order_>=0)
      { ifm >>= order_; ifp >>= order_; }
    else
      { ifm /= nside_; ifp /= nside_; }
    face_num = (ifp==ifm) ? ((ifp&3)+4) : ((ifp<ifm) ? ifp : (ifm+8));
    }
  else // South Polar cap
    {
    int ip = npix_ - pix;
    iring = (1+isqrt(2*ip-1))>>1; //counted from South pole
    iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));
    kshift = 0;
    nr = iring;
    iring = 2*nl2-iring;
    face_num = 8 + (iphi-1)/nr;
    }

  int irt = iring - (jrll[face_num]*nside_) + 1;
  int ipt = 2*iphi- jpll[face_num]*nr - kshift -1;
  if (ipt>=nl2) ipt-=8*nside_;

  ix =  (ipt-irt) >>1;
  iy = (-ipt-irt) >>1;
  }

int Healpix_Base::xyf2ring (int ix, int iy, int face_num) const
  {
  int nl4 = 4*nside_;
  int jr = (jrll[face_num]*nside_) - ix - iy  - 1;

  int nr, kshift, n_before;
  if (jr<nside_)
    {
    nr = jr;
    n_before = 2*nr*(nr-1);
    kshift = 0;
    }
  else if (jr > 3*nside_)
    {
    nr = nl4-jr;
    n_before = npix_ - 2*(nr+1)*nr;
    kshift = 0;
    }
  else
    {
    nr = nside_;
    n_before = ncap_ + (jr-nside_)*nl4;
    kshift = (jr-nside_)&1;
    }

  int jp = (jpll[face_num]*nr + ix - iy + 1 + kshift) / 2;
  planck_assert(jp<=4*nr,"must not happen");
  if (jp<1) jp+=nl4; // assumption: if this triggers, then nl4==4*nr

  return n_before + jp - 1;
  }

Healpix_Base::Healpix_Base ()
  : order_(-1), nside_(0), npface_(0), ncap_(0), npix_(0),
    fact1_(0), fact2_(0), scheme_(RING) {}

void Healpix_Base::Set (int order, Healpix_Ordering_Scheme scheme)
  {
  planck_assert ((order>=0)&&(order<=order_max), "bad order");
  order_  = order;
  nside_  = 1<<order;
  npface_ = nside_<<order_;
  ncap_   = (npface_-nside_)<<1;
  npix_   = 12*npface_;
  fact2_  = 4./npix_;
  fact1_  = (nside_<<1)*fact2_;
  scheme_ = scheme;
  }
void Healpix_Base::SetNside (int nside, Healpix_Ordering_Scheme scheme)
  {
  order_  = nside2order(nside);
  planck_assert ((scheme!=NEST) || (order_>=0),
    "SetNside: nside must be power of 2 for nested maps");
  nside_  = nside;
  npface_ = nside_*nside_;
  ncap_   = (npface_-nside_)<<1;
  npix_   = 12*npface_;
  fact2_  = 4./npix_;
  fact1_  = (nside_<<1)*fact2_;
  scheme_ = scheme;
  }

double Healpix_Base::ring2z (int ring) const
  {
  if (ring<nside_)
    return 1 - ring*ring*fact2_;
  if (ring <=3*nside_)
    return (2*nside_-ring)*fact1_;
  ring=4*nside_ - ring;
  return ring*ring*fact2_ - 1;
  }

int Healpix_Base::pix2ring (int pix) const
  {
  if (scheme_==RING)
    {
    if (pix<ncap_) // North Polar cap
      return (1+isqrt(1+2*pix))>>1; //counted from North pole
    else if (pix<(npix_-ncap_)) // Equatorial region
      return (pix-ncap_)/(4*nside_) + nside_; // counted from North pole
    else // South Polar cap
      return 4*nside_-((1+isqrt(2*(npix_-pix)-1))>>1); //counted from South pole
    }
  else
    {
    int face_num, ix, iy;
    nest2xyf(pix,ix,iy,face_num);
    return (jrll[face_num]<<order_) - ix - iy - 1;
    }
  }

int Healpix_Base::nest2ring (int pix) const
  {
  planck_assert(order_>=0, "nest2ring: need hierarchical map");
  int ix, iy, face_num;
  nest2xyf (pix, ix, iy, face_num);
  return xyf2ring (ix, iy, face_num);
  }

int Healpix_Base::ring2nest (int pix) const
  {
  planck_assert(order_>=0, "ring2nest: need hierarchical map");
  int ix, iy, face_num;
  ring2xyf (pix, ix, iy, face_num);
  return xyf2nest (ix, iy, face_num);
  }

int Healpix_Base::nest_peano_helper (int pix, int dir) const
  {
  int face = pix>>(2*order_), result = 0;
  uint8 path = peano_face2path[dir][face];

  for (int shift=2*order_-2; shift>=0; shift-=2)
    {
    uint8 spix = uint8((pix>>shift) & 0x3);
    result <<= 2;
    result |= peano_subpix[dir][path][spix];
    path=peano_subpath[dir][path][spix];
    }

  return result + (int(peano_face2face[dir][face])<<(2*order_));
  }

int Healpix_Base::nest2peano (int pix) const
  { return nest_peano_helper(pix,0); }

int Healpix_Base::peano2nest (int pix) const
  { return nest_peano_helper(pix,1); }

int Healpix_Base::zphi2pix (double z, double phi) const
  {
  double za = abs(z);
  double tt = fmodulo(phi*inv_halfpi,4.0); // in [0,4)

  if (scheme_==RING)
    {
    if (za<=twothird) // Equatorial region
      {
      unsigned int nl4 = 4*nside_;
      double temp1 = nside_*(0.5+tt);
      double temp2 = nside_*z*0.75;
      int jp = int(temp1-temp2); // index of  ascending edge line
      int jm = int(temp1+temp2); // index of descending edge line

      // ring number counted from z=2/3
      int ir = nside_ + 1 + jp - jm; // in {1,2n+1}
      int kshift = 1-(ir&1); // kshift=1 if ir even, 0 otherwise

      unsigned int t1 = jp+jm-nside_+kshift+1+nl4+nl4;
      int ip = (order_>0) ?
        (t1>>1)&(nl4-1) : ((t1>>1)%nl4); // in {0,4n-1}

      return ncap_ + (ir-1)*nl4 + ip;
      }
    else  // North & South polar caps
      {
      double tp = tt-int(tt);
      double tmp = nside_*sqrt(3*(1-za));

      int jp = int(tp*tmp); // increasing edge line index
      int jm = int((1.0-tp)*tmp); // decreasing edge line index

      int ir = jp+jm+1; // ring number counted from the closest pole
      int ip = int(tt*ir); // in {0,4*ir-1}
      planck_assert((ip>=0) &&(ip<4*ir),"must not happen");
      //ip %= 4*ir;

      return (z>0)  ?  2*ir*(ir-1) + ip  :  npix_ - 2*ir*(ir+1) + ip;
      }
    }
  else // scheme_ == NEST
    {
    if (za<=twothird) // Equatorial region
      {
      double temp1 = nside_*(0.5+tt);
      double temp2 = nside_*(z*0.75);
      int jp = int(temp1-temp2); // index of  ascending edge line
      int jm = int(temp1+temp2); // index of descending edge line
      int ifp = jp >> order_;  // in {0,4}
      int ifm = jm >> order_;
      int face_num = (ifp==ifm) ? ((ifp&3)+4) : ((ifp<ifm) ? ifp : (ifm+8));

      int ix = jm & (nside_-1),
          iy = nside_ - (jp & (nside_-1)) - 1;
      return xyf2nest(ix,iy,face_num);
      }
    else // polar region, za > 2/3
      {
      int ntt = min(3,int(tt));
      double tp = tt-ntt;
      double tmp = nside_*sqrt(3*(1-za));

      int jp = int(tp*tmp); // increasing edge line index
      int jm = int((1.0-tp)*tmp); // decreasing edge line index
      if (jp>=nside_) jp = nside_-1; // for points too close to the boundary
      if (jm>=nside_) jm = nside_-1;
      return (z>=0) ?
        xyf2nest(nside_-jm -1,nside_-jp-1,ntt) :
        xyf2nest(jp,jm,ntt+8);
      }
    }
  }

void Healpix_Base::pix2zphi (int pix, double &z, double &phi) const
  {
  if (scheme_==RING)
    {
    if (pix<ncap_) // North Polar cap
      {
      int iring = (1+isqrt(1+2*pix))>>1; //counted from North pole
      int iphi  = (pix+1) - 2*iring*(iring-1);

      z = 1.0 - (iring*iring)*fact2_;
      phi = (iphi-0.5) * halfpi/iring;
      }
    else if (pix<(npix_-ncap_)) // Equatorial region
      {
      int nl4 = 4*nside_;
      int ip  = pix - ncap_;
      int tmp = (order_>=0) ? ip>>(order_+2) : ip/nl4;
      int iring = tmp + nside_,
          iphi = ip-nl4*tmp+1;;
      // 1 if iring+nside is odd, 1/2 otherwise
      double fodd = ((iring+nside_)&1) ? 1 : 0.5;

      z = (2*nside_-iring)*fact1_;
      phi = (iphi-fodd) * pi*0.75*fact1_;
      }
    else // South Polar cap
      {
      int ip = npix_ - pix;
      int iring = (1+isqrt(2*ip-1))>>1; //counted from South pole
      int iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1));

      z = -1.0 + (iring*iring)*fact2_;
      phi = (iphi-0.5) * halfpi/iring;
      }
    }
  else
    {
    int face_num, ix, iy;
    nest2xyf(pix,ix,iy,face_num);

    int jr = (jrll[face_num]<<order_) - ix - iy - 1;

    int nr;
    if (jr<nside_)
      {
      nr = jr;
      z = 1 - nr*nr*fact2_;
      }
    else if (jr > 3*nside_)
      {
      nr = nside_*4-jr;
      z = nr*nr*fact2_ - 1;
      }
    else
      {
      nr = nside_;
      z = (2*nside_-jr)*fact1_;
      }

    int tmp=jpll[face_num]*nr+ix-iy;
    planck_assert(tmp<8*nr,"must not happen");
    if (tmp<0) tmp+=8*nr;
    phi = (nr==nside_) ? 0.75*halfpi*tmp*fact1_ :
                         (0.5*halfpi*tmp)/nr;
    }
  }

void Healpix_Base::query_polygon (const vector<pointing> &vertex,
  bool inclusive, rangeset<int> &pixset) const
  {
  tsize nv=vertex.size();
  tsize ncirc = inclusive ? nv+1 : nv;
  planck_assert(nv>=3,"not enough vertices in polygon");
  arr<vec3> vv(nv);
  for (tsize i=0; i<nv; ++i)
    vv[i]=vertex[i].to_vec3();
  arr<vec3> normal(ncirc);
  int flip=0;
  for (tsize i=0; i<nv; ++i)
    {
    normal[i]=crossprod(vv[i],vv[(i+1)%nv]);
    double hnd=dotprod(normal[i],vv[(i+2)%nv]);
    planck_assert(abs(hnd)>1e-10,"degenerate corner");
    if (i==0)
      flip = (hnd<0.) ? -1 : 1;
    else
      planck_assert(flip*hnd>0,"polygon is not convex");
    normal[i]*=flip/normal[i].Length();
    }
  arr<double> rad(ncirc,halfpi);
  if (inclusive)
    {
    double cosrad;
    find_enclosing_circle (vv, normal[nv], cosrad);
    rad[nv]=acos(cosrad);
    }
  query_multidisc(normal,rad,inclusive,pixset);
  }

void Healpix_Base::get_ring_info_small (int ring, int &startpix, int &ringpix,
  bool &shifted) const
  {
  int northring = (ring>2*nside_) ? 4*nside_-ring : ring;
  if (northring < nside_)
    {
    shifted = true;
    ringpix = 4*northring;
    startpix = 2*northring*(northring-1);
    }
  else
    {
    shifted = ((northring-nside_) & 1) == 0;
    ringpix = 4*nside_;
    startpix = ncap_ + (northring-nside_)*ringpix;
    }
  if (northring != ring) // southern hemisphere
    startpix = npix_ - startpix - ringpix;
  }
void Healpix_Base::get_ring_info (int ring, int &startpix, int &ringpix,
  double &costheta, double &sintheta, bool &shifted) const
  {
  int northring = (ring>2*nside_) ? 4*nside_-ring : ring;
  if (northring < nside_)
    {
    double tmp = northring*northring*fact2_;
    costheta = 1 - tmp;
    sintheta = sqrt(tmp*(2-tmp));
    ringpix = 4*northring;
    shifted = true;
    startpix = 2*northring*(northring-1);
    }
  else
    {
    costheta = (2*nside_-northring)*fact1_;
    sintheta = sqrt((1+costheta)*(1-costheta));
    ringpix = 4*nside_;
    shifted = ((northring-nside_) & 1) == 0;
    startpix = ncap_ + (northring-nside_)*ringpix;
    }
  if (northring != ring) // southern hemisphere
    {
    costheta = -costheta;
    startpix = npix_ - startpix - ringpix;
    }
  }
void Healpix_Base::get_ring_info2 (int ring, int &startpix, int &ringpix,
  double &theta, bool &shifted) const
  {
  int northring = (ring>2*nside_) ? 4*nside_-ring : ring;
  if (northring < nside_)
    {
    double tmp = northring*northring*fact2_;
    double costheta = 1 - tmp;
    double sintheta = sqrt(tmp*(2-tmp));
    theta = atan2(sintheta,costheta);
    ringpix = 4*northring;
    shifted = true;
    startpix = 2*northring*(northring-1);
    }
  else
    {
    theta = acos((2*nside_-northring)*fact1_);
    ringpix = 4*nside_;
    shifted = ((northring-nside_) & 1) == 0;
    startpix = ncap_ + (northring-nside_)*ringpix;
    }
  if (northring != ring) // southern hemisphere
    {
    theta = pi-theta;
    startpix = npix_ - startpix - ringpix;
    }
  }

void Healpix_Base::neighbors (int pix, fix_arr<int,8> &result) const
  {
  int ix, iy, face_num;
  (scheme_==RING) ?
    ring2xyf(pix,ix,iy,face_num) : nest2xyf(pix,ix,iy,face_num);

  const int nsm1 = nside_-1;
  if ((ix>0)&&(ix<nsm1)&&(iy>0)&&(iy<nsm1))
    {
    if (scheme_==RING)
      for (int m=0; m<8; ++m)
        result[m] = xyf2ring(ix+nb_xoffset[m],iy+nb_yoffset[m],face_num);
    else
      {
      int fpix = face_num<<(2*order_),
          px0=spread_bits(ix  ), py0=spread_bits(iy  )<<1,
          pxp=spread_bits(ix+1), pyp=spread_bits(iy+1)<<1,
          pxm=spread_bits(ix-1), pym=spread_bits(iy-1)<<1;

      result[0] = fpix+pxm+py0; result[1] = fpix+pxm+pyp;
      result[2] = fpix+px0+pyp; result[3] = fpix+pxp+pyp;
      result[4] = fpix+pxp+py0; result[5] = fpix+pxp+pym;
      result[6] = fpix+px0+pym; result[7] = fpix+pxm+pym;
      }
    }
  else
    {
    for (int i=0; i<8; ++i)
      {
      int x=ix+nb_xoffset[i], y=iy+nb_yoffset[i];
      int nbnum=4;
      if (x<0)
        { x+=nside_; nbnum-=1; }
      else if (x>=nside_)
        { x-=nside_; nbnum+=1; }
      if (y<0)
        { y+=nside_; nbnum-=3; }
      else if (y>=nside_)
        { y-=nside_; nbnum+=3; }

      int f = nb_facearray[nbnum][face_num];
      if (f>=0)
        {
        int bits = nb_swaparray[nbnum][face_num>>2];
        if (bits&1) x=nside_-x-1;
        if (bits&2) y=nside_-y-1;
        if (bits&4) std::swap(x,y);
        result[i] = (scheme_==RING) ? xyf2ring(x,y,f) : xyf2nest(x,y,f);
        }
      else
        result[i] = -1;
      }
    }
  }

void Healpix_Base::get_interpol (const pointing &ptg, fix_arr<int,4> &pix,
  fix_arr<double,4> &wgt) const
  {
  double z = cos (ptg.theta);
  int ir1 = ring_above(z);
  int ir2 = ir1+1;
  double theta1, theta2, w1, tmp, dphi;
  int sp,nr;
  bool shift;
  int i1,i2;
  if (ir1>0)
    {
    get_ring_info2 (ir1, sp, nr, theta1, shift);
    dphi = twopi/nr;
    tmp = (ptg.phi/dphi - .5*shift);
    i1 = (tmp<0) ? int(tmp)-1 : int(tmp);
    w1 = (ptg.phi-(i1+.5*shift)*dphi)/dphi;
    i2 = i1+1;
    if (i1<0) i1 +=nr;
    if (i2>=nr) i2 -=nr;
    pix[0] = sp+i1; pix[1] = sp+i2;
    wgt[0] = 1-w1; wgt[1] = w1;
    }
  if (ir2<(4*nside_))
    {
    get_ring_info2 (ir2, sp, nr, theta2, shift);
    dphi = twopi/nr;
    tmp = (ptg.phi/dphi - .5*shift);
    i1 = (tmp<0) ? int(tmp)-1 : int(tmp);
    w1 = (ptg.phi-(i1+.5*shift)*dphi)/dphi;
    i2 = i1+1;
    if (i1<0) i1 +=nr;
    if (i2>=nr) i2 -=nr;
    pix[2] = sp+i1; pix[3] = sp+i2;
    wgt[2] = 1-w1; wgt[3] = w1;
    }

  if (ir1==0)
    {
    double wtheta = ptg.theta/theta2;
    wgt[2] *= wtheta; wgt[3] *= wtheta;
    double fac = (1-wtheta)*0.25;
    wgt[0] = fac; wgt[1] = fac; wgt[2] += fac; wgt[3] +=fac;
    pix[0] = (pix[2]+2)%4;
    pix[1] = (pix[3]+2)%4;
    }
  else if (ir2==4*nside_)
    {
    double wtheta = (ptg.theta-theta1)/(pi-theta1);
    wgt[0] *= (1-wtheta); wgt[1] *= (1-wtheta);
    double fac = wtheta*0.25;
    wgt[0] += fac; wgt[1] += fac; wgt[2] = fac; wgt[3] =fac;
    pix[2] = ((pix[0]+2)&3)+npix_-4;
    pix[3] = ((pix[1]+2)&3)+npix_-4;
    }
  else
    {
    double wtheta = (ptg.theta-theta1)/(theta2-theta1);
    wgt[0] *= (1-wtheta); wgt[1] *= (1-wtheta);
    wgt[2] *= wtheta; wgt[3] *= wtheta;
    }

  if (scheme_==NEST)
    for (tsize m=0; m<pix.size(); ++m)
      pix[m] = ring2nest(pix[m]);
  }

void Healpix_Base::swap (Healpix_Base &other)
  {
  std::swap(order_,other.order_);
  std::swap(nside_,other.nside_);
  std::swap(npface_,other.npface_);
  std::swap(ncap_,other.ncap_);
  std::swap(npix_,other.npix_);
  std::swap(fact1_,other.fact1_);
  std::swap(fact2_,other.fact2_);
  std::swap(scheme_,other.scheme_);
  }

double Healpix_Base::max_pixrad() const
  {
  vec3 va,vb;
  va.set_z_phi (2./3., pi/(4*nside_));
  double t1 = 1.-1./nside_;
  t1*=t1;
  vb.set_z_phi (1-t1/3, 0);
  return v_angle(va,vb);
  }

double Healpix_Base::max_pixrad(int ring) const
  {
  if (ring>=2*nside_) ring=4*nside_-ring;
  double z=ring2z(ring), z_up=(ring>1) ? ring2z(ring-1) : 1.;
  vec3 mypos, uppos;
  uppos.set_z_phi(z_up,0);
  if (ring<=nside_)
    {
    mypos.set_z_phi(z,pi/(4*ring));
    return v_angle(mypos,uppos);
    }
  mypos.set_z_phi(z,0);
  double vdist=v_angle(mypos,uppos);
  double hdist=sqrt(1.-z*z)*pi/(4*nside_);
  return max(hdist,vdist);
  }

arr<int> Healpix_Base::swap_cycles() const
  {
  planck_assert(order_>=0, "swap_cycles(): need hierarchical map");
  arr<int> result(swap_clen[order_]);
  tsize ofs=0;
  for (int m=0; m<order_;++m) ofs+=swap_clen[m];
  for (tsize m=0; m<result.size();++m) result[m]=swap_cycle[m+ofs];
  return result;
  }
