#ifndef PBC_H
#define PBC_H

#include "Vector3D.hpp"

std::vector<Vector3D> shifts_for_pbc;

void pbc_setup(const Vector3D& cryst)
{
   shifts_for_pbc.clear();
   shifts_for_pbc.reserve(27);
   for(int x=-1;x<=+1;++x)
   for(int y=-1;y<=+1;++y)
   for(int z=-1;z<=+1;++z)
   {
      shifts_for_pbc.push_back(Vector3D(x*cryst.x,y*cryst.y,z*cryst.z));
   }
   return ;
}

Vector3D pbc(Vector3D v, const Vector3D& ref)
{
   double min=DBL_MAX;
   Vector3D res;
   Vector3D tmp;
   for(size_t i=0,i_size=shifts_for_pbc.size();i<i_size;++i)
   {
      tmp=v+shifts_for_pbc.at(i);
      const double dis=(tmp-ref).norm2();
      if(min>dis)
      {
         min=dis;
         res=tmp;   
      }
   }
   return res;
}

#endif 

