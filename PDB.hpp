#ifndef PDB_H
#define PDB_H

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <random>
#include "Vector3D.hpp"
#include "ReadFile.hpp"
#include "PBC.hpp"
#include "Kahan_summation_algorithm.hpp"
#include "Quaternion.hpp"

std::vector<std::string> split(const std::string req);

class Atom;
class Protein;

std::mt19937_64 mt(472354627351723);

class Protein
{
   public:
   Vector3D cryst;
   double time;
   std::vector<Atom> atoms;   
   public:
   Protein(){}
   public:
   int get_index_of_first_resudue()const;
   int get_index_of_last_resudue()const;
   std::vector<Vector3D> get_atoms_of_residue(const int& index)const;
   Vector3D get_center_of_residue(const int& index)const;
   std::vector<Vector3D> get_CAs()const;
   void fit_to_(const Protein& ref, const int& index_res, const int step_trial=100000);
};

class Atom
{
   public:
   int index_atom;       //index of atom
   std::string pt;       //type of particle
   std::string res;      //name of residue
   std::string domain;   //name of domain
   int index_res;        //index of residue
   Vector3D position;     //position of atom
   double occupancy;
   double b_factor;
   Atom(const std::vector<std::string>& vs);
};


Atom::Atom(const std::vector<std::string>& vs)
{
   index_atom = boost::lexical_cast<int> (vs.at(1));
   pt         = vs.at(2);
   res        = vs.at(3);
   index_res  = boost::lexical_cast<int> (vs.at(4));
   position    = Vector3D(boost::lexical_cast<double> (vs.at(5)),boost::lexical_cast<double> (vs.at(6)),boost::lexical_cast<double> (vs.at(7)));
   occupancy  = boost::lexical_cast<double> (vs.at(8));
   b_factor   = boost::lexical_cast<double> (vs.at(9));
}

int Protein::get_index_of_first_resudue()const
{
     return (atoms.front()).index_res; 
}

int Protein::get_index_of_last_resudue()const
{
   return (atoms.back()).index_res; 
}

std::vector<Vector3D> Protein::get_atoms_of_residue(const int& index)const
{
   std::vector<Vector3D> result;
   for(size_t i=0,i_size=atoms.size();i<i_size;++i)
   {
      const Atom& a = atoms.at(i);
      if(index==a.index_res && "CA"!=a.pt)
      {
         result.push_back(a.position);
      }
   }
   return result; 
}

Vector3D Protein::get_center_of_residue(const int& index)const
{
   std::vector<Vector3D> ps = get_atoms_of_residue(index);
   for(size_t i=1,i_size=ps.size();i<i_size;++i)
   {
      ps.at(i)=pbc(ps.at(i),ps.at(0));   
   }
   Vector3D result;
   for(size_t i=0,i_size=ps.size();i<i_size;++i)
   {
      result += ps.at(i);   
   }
   result/=ps.size();
   return result;
}

std::vector<Vector3D> Protein::get_CAs()const
{
   std::vector<Vector3D> res;
   for(size_t i=0,i_size=atoms.size();i<i_size;++i)
   {
      if("CA"==atoms.at(i).pt) 
      {
         res.push_back(atoms.at(i).position);
      }
   }
   return res;
}

void Protein::fit_to_(const Protein& ref, const int& index_res, const int step_trial)
{

   const auto centering =
      [](std::vector<Vector3D> res, const Vector3D& center)->std::tuple<Vector3D,std::vector<Vector3D> >
      {
         //Kahan x;
         //Kahan y;
         //Kahan z;
         //for(size_t i=0,i_size=res.size();i<i_size;++i)
         //{
         //   x += res.at(i).x;
         //   y += res.at(i).y;
         //   z += res.at(i).z;
         //}
         //const Vector3D center(x.get_av(),y.get_av(),z.get_av());
         for(size_t i=0,i_size=res.size();i<i_size;++i)
         {
            res.at(i) -= center;
         }
         return std::tuple<Vector3D,std::vector<Vector3D> >(center,res);
      };
   const auto refca  = centering(ref.get_CAs(),ref.get_center_of_residue(index_res));
   const auto rmsd_CA = [&](const auto& tmpca)
   {
      Kahan res=0.0;
      for(size_t i=0,i_size=(std::get<1>(refca)).size();i<i_size;++i)
      {
         res+=(tmpca.at(i)-(std::get<1>(refca)).at(i)).norm2();
      }
      return res.get_av();
   };

   Vector3D rot_current;
   Vector3D rot_tmp;
   std::tuple<Vector3D,std::vector<Vector3D> > struct_begin  = centering(get_CAs(),get_center_of_residue(index_res));
   std::vector<Vector3D> struct_current=std::get<1>(struct_begin);
   std::vector<Vector3D> struct_tmp    =std::get<1>(struct_begin);
   double rmsd_current=rmsd_CA(std::get<1>(struct_begin));
   double rmsd_tmp;
   constexpr double DELTA_RANGE = M_PI/50.0;
   const Vector3D ax(1.0,0.0,0.0);
   const Vector3D ay(0.0,1.0,0.0);
   const Vector3D az(0.0,0.0,1.0);
   for(int step=0;step<step_trial;++step)
   {
      Vector3D delta;
      std::uniform_real_distribution<double> dist(-DELTA_RANGE,+DELTA_RANGE); 
      do
      {
         delta.x = dist(mt);
         delta.y = dist(mt);
         delta.z = dist(mt);
      }while(delta.norm2()>DELTA_RANGE);
      rot_tmp=rot_current+delta;
      struct_tmp = [&]()
         {
            std::vector<Vector3D> res(struct_current.size());
            Quaternion qx(ax,rot_tmp.x);
            Quaternion qy(ay,rot_tmp.y);
            Quaternion qz(az,rot_tmp.z);
            for(size_t i=0,i_size=res.size();i<i_size;++i)
            {
               Quaternion q = qz*qy*qx*(std::get<1>(struct_begin)).at(i)*bar(qx)*bar(qy)*bar(qz);
               res.at(i) = Vector3D(q.b, q.c, q.d);
            }
            return res;
         }();
      rmsd_tmp = rmsd_CA(struct_tmp);
      if(rmsd_current>rmsd_tmp)
      {
         rot_current    = rot_tmp;
         struct_current = struct_tmp;
         rmsd_current   = rmsd_tmp;
      }
   }
   {//apply the result
      Quaternion qx(ax,rot_current.x);
      Quaternion qy(ay,rot_current.y);
      Quaternion qz(az,rot_current.z);
      for(size_t i=0,i_size=atoms.size();i<i_size;++i)
      {
         Vector3D& p = atoms.at(i).position;
         p -= std::get<0>(struct_begin);
         Quaternion q = qz*qy*qx*p*bar(qx)*bar(qy)*bar(qz);
         p = Vector3D(q.b, q.c, q.d);
         p += std::get<0>(refca);
      }
   }
}

std::vector<Protein> load(const std::string file_name, const std::tuple<double,double>& t_range)
{
   std::vector<Protein> result;
   Getline gl(file_name);
   std::string tmp;
   Protein ptmp;
   while(gl.is_open())
   {
      try{
      tmp = gl.get(); 
      }catch(...){break;}
      std::vector<std::string> vs = split(tmp);
      if("TITLE"==vs.at(0))
      {
         ptmp.time=boost::lexical_cast<double>(vs.back());
      }
      if("CRYST1"==vs.at(0))
      {
         ptmp.cryst = Vector3D(boost::lexical_cast<double>(vs.at(1)),boost::lexical_cast<double>(vs.at(2)),boost::lexical_cast<double>(vs.at(3)));
      }
      if("ATOM"==vs.at(0))
      {
         ptmp.atoms.push_back(Atom(vs));
      }
      if("ENDMDL"==vs.at(0))
      {
         if(std::get<0>(t_range)<=ptmp.time && ptmp.time<=std::get<1>(t_range))
         {
            result.push_back(ptmp);
         }
      }
   }
   return result;
}

Protein get_next(Getline& gl, const std::tuple<double,double>& t_range, bool is_first=false)
{
   Protein res;
   std::string tmp;
   while(gl.is_open())
   {
      try{
      tmp = gl.get(); 
      }catch(...){break;}
      std::vector<std::string> vs = split(tmp);
      if("TITLE"==vs.at(0))
      {
         res.time=boost::lexical_cast<double>(vs.back());
      }
      if("CRYST1"==vs.at(0))
      {
         res.cryst = Vector3D(boost::lexical_cast<double>(vs.at(1)),boost::lexical_cast<double>(vs.at(2)),boost::lexical_cast<double>(vs.at(3)));
      }
      if("ATOM"==vs.at(0))
      {
         res.atoms.push_back(Atom(vs));
      }
      if("ENDMDL"==vs.at(0))
      {
         if(std::get<0>(t_range)<=res.time && res.time<=std::get<1>(t_range))
         {
            return res;
         }
         else
         {
            if(!is_first)
            {
               throw "over_run";
            }
         }
      }
   }
   throw "error";
   return Protein();
}

std::vector<std::string> split(const std::string req)
{
   std::vector<std::string> tmp_vsplit;
   boost::algorithm::split(tmp_vsplit,req,boost::is_any_of(" "));
   std::vector<std::string> vsplit;
   for(int i=0,size=tmp_vsplit.size();i<size;++i)
   {
      if(""!=tmp_vsplit[i])
      {
         vsplit.push_back(tmp_vsplit[i]);
      }
   }
   return vsplit;
}

#endif
