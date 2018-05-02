#ifndef PDB_H
#define PDB_H

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "Vector3D.hpp"
#include "ReadFile.hpp"

std::vector<std::string> split(const std::string req);

class Atom;
class Protein;

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
   domain     = vs.at(4);
   index_res  = boost::lexical_cast<int> (vs.at(5));
   position    = Vector3D(boost::lexical_cast<int> (vs.at(6)),boost::lexical_cast<int> (vs.at(7)),boost::lexical_cast<int> (vs.at(8)));
   occupancy  = boost::lexical_cast<int> (vs.at(9));
   b_factor   = boost::lexical_cast<int> (vs.at(10));
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
   const std::vector<Vector3D> ps = get_atoms_of_residue(index);
   Vector3D result;
   for(size_t i=0,i_size=ps.size();i<i_size;++i)
   {
      result += ps.at(i);   
   }
   result/=ps.size();
   return result;
}

std::vector<Protein> load(const std::string file_name)
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
         result.push_back(ptmp);
      }
   }
   return result;
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
