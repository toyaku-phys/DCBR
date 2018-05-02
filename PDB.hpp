#ifndef PDB_H
#define PDB_H

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "Vector3D.hpp"
#include "ReadFile.hpp"
//ENDMDLできれるお

std::vector<std::string> split(const std::string req);

class Atom;

class Protein
{
   public:
   Vector3D cryst;
   double time;
   std::vector<Atom> atoms;   
   public:
   Protein(){}
};

class Atom
{
   public:
   int index_atom;       //index of atom
   std::string pt;       //type of particle
   std::string res;      //name of residue
   std::string domain;   //name of domain
   int index_res;        //index of residue
   Vector3D positon;     //positon of atom
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
   positon    = Vector3D(boost::lexical_cast<int> (vs.at(6)),boost::lexical_cast<int> (vs.at(7)),boost::lexical_cast<int> (vs.at(8)));
   occupancy  = boost::lexical_cast<int> (vs.at(9));
   b_factor   = boost::lexical_cast<int> (vs.at(10));
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
