#include "PDB.hpp"
#include "Kahan_summation_algorithm.hpp"

Vector3D displacement_of_residue(const Protein& protein_a,const Protein& protein_b,const int& index_target_residue);
Vector3D velocity_of_residue
(const Protein& protein_a,const Protein& protein_b,const int& index_target_residue);
double correlation(const Vector3D& a, const Vector3D& b);
void plot_ca_for_yuba(const std::string file_name, const int step, const Protein& pr);

int main(int argc, char* argv[])
{
   std::vector<std::string> argments;
   for(int i=0;i<argc;++i){ argments.push_back(argv[i]); }

   std::string input_file_name;
   std::string output_file_name;
   int target_residue_index;
   std::tuple<double,double> time_range(-DBL_MAX,+DBL_MAX);
   bool silent=false;
   std::tuple<int,int> delta_range;

   std::bitset<5> flags("00000");

   try
   {
      for(size_t i=0,i_size=argments.size();i<i_size;++i)
      {
         std::vector<std::string> vs;
         boost::algorithm::split(vs,argments.at(i),boost::is_any_of("="));
         if("in"==vs.at(0))
         {
            input_file_name=vs.at(1);
            flags.set(0);
         }
         if("out"==vs.at(0))
         {
            output_file_name=vs.at(1);
            flags.set(1);
         }
         if("target"==vs.at(0))
         {
            target_residue_index=boost::lexical_cast<int>(vs.at(1));
            flags.set(2);
         }
         if("time"==vs.at(0))
         {
            std::vector<std::string> vs_;
            boost::algorithm::split(vs_,vs.at(1),boost::is_any_of("-"));
            time_range=std::tuple<double,double>
            (
               boost::lexical_cast<double>(vs_.at(0)),
               boost::lexical_cast<double>(vs_.at(1))
            ); 
            flags.set(3);
         }
         if("delta_range"==vs.at(0))
         {
            std::vector<std::string> vs_;
            boost::algorithm::split(vs_,vs.at(1),boost::is_any_of("-"));
            std::get<0>(delta_range) = boost::lexical_cast<int> (vs_.at(0));
            std::get<1>(delta_range) = boost::lexical_cast<int> (vs_.at(1));
            flags.set(4);
         }
         if("--silent"==argments.at(i))
         {
            silent=true;
         }
      }
      std::cout<<flags.count()<<std::endl<<std::flush; 
      if(!flags.all()){throw "flag";}
   }catch(...)
   {
      std::cout<<"Failure @input argments"<<std::endl;exit(0);
   }
   
   int rbegn = 1;
   std::vector<std::vector<Kahan> > correlationss;
   bool f_delta_t = false;
   double Delta_T = 0.0;
   for(int d=std::get<0>(delta_range);d<=std::get<1>(delta_range);++d)
   {
      Getline gl_m(input_file_name);
      Getline gl_z(input_file_name);
      std::tuple<Protein,Protein> minus(get_next(gl_m,time_range,true),get_next(gl_m,time_range));//-+
      std::tuple<Protein,Protein> zero(get_next(gl_z,time_range,true),get_next(gl_z,time_range));
      //shift delta
      for(int dlt=0;dlt<d;++dlt)
      {
         Protein tmp = std::get<1> (zero);
         std::get<1> (zero) = get_next(gl_z,time_range);
         if(!f_delta_t){Delta_T=std::abs(std::get<1>(zero).time-tmp.time);}
         std::get<0> (zero) = tmp;
      }
                rbegn=(std::get<0>(zero)).get_index_of_first_resudue();
      const int rlast=(std::get<0>(zero)).get_index_of_last_resudue();
      std::vector<Kahan> correlations(rlast-rbegn+1);
      for(;;)
      {
         //each res
         const Vector3D mv_pos = velocity_of_residue(std::get<0>(minus),std::get<1>(minus),target_residue_index);
         for(int r=rbegn;r<=rlast;++r)
         {
            const Vector3D v_j = velocity_of_residue(std::get<0>(zero),std::get<1>(zero),r);
            correlations.at(r-rbegn) += correlation(mv_pos,v_j)/(mv_pos.norm()*v_j.norm());
         }
         try
         {
            Protein tmp = std::get<1>(zero); 
            std::get<1> (zero) = get_next(gl_z,time_range);
            std::get<0> (zero) = tmp;
            tmp = std::get<1>(minus);
            std::get<1> (minus) = get_next(gl_m,time_range);
            std::get<0> (minus) = tmp;
         }catch(...){break;}
         if(!silent)
         {
            static int idx=0;
            switch(idx++)
            {
               case 0: std::cout<<"ま"<<std::flush;break;
               case 1: 
               case 2: std::cout<<"だ"<<std::flush;break;
               case 3: std::cout<<"よ"<<std::flush;break;
               default: std::cout<<"!"<<std::flush;idx=0;break;
            };
         }
      }//end of ;;
      correlationss.push_back(correlations);
   }

   {
      std::ofstream ofs(output_file_name,std::ios::trunc);
      int idx=0;
      ofs<<"# n C(R_1) C(R_2) ..."<<std::endl;
      for(int d=std::get<0>(delta_range);d<=std::get<1>(delta_range);++d,++idx)
      {
         ofs<<d*Delta_T<<" ";
         const std::vector<Kahan>& rs = correlationss.at(idx);
         for(size_t i=0,i_size=rs.size();i<i_size;++i)
         {
            ofs<<rs.at(i).get_av()<<" ";
         }
      }
      
   }

   if(!silent)
   {
      std::cout<<std::endl;
      std::cout<<"もういいよ！"<<std::endl;
   }

   return EXIT_SUCCESS;
}

Vector3D displacement_of_residue
(
   const Protein& protein_a,
   const Protein& protein_b,
   const int& index_target_residue
)
{
   const Vector3D b = protein_b.get_center_of_residue(index_target_residue);
   const Vector3D a = pbc(protein_a.get_center_of_residue(index_target_residue),b);
   return b-a;
}

Vector3D velocity_of_residue
(
   const Protein& protein_a,
   const Protein& protein_b,
   const int& index_target_residue
)
{
   const Vector3D b = protein_b.get_center_of_residue(index_target_residue);
   const Vector3D a = pbc(protein_a.get_center_of_residue(index_target_residue),b);
   return (b-a)/(protein_b.time-protein_a.time);
}


double correlation(const Vector3D& u_i, const Vector3D& u_j)
{
   return u_i*u_j;   
}

void plot_ca_for_yuba(const std::string file_name, const int step, const Protein& pr)
{
   const std::vector<Vector3D> cas = pr.get_CAs();
   if(0==step)
   {
      std::ofstream ofs(file_name,std::ios::trunc);
      ofs<<"# Coordinate "<<step<<std::endl;
      for(size_t i=0,i_size=cas.size();i<i_size;++i)
      {
         ofs<<i<<" "<<cas.at(i)<<std::endl;
      }
      ofs<<std::endl;
      ofs<<"# Bond "<<step<<std::endl;
      for(size_t i=1,i_size=cas.size();i<i_size;++i)
      {
         ofs<<i-1<<" "<<i<<std::endl;
      }
      ofs<<std::endl<<std::endl;
      ofs.close();
   }
   else
   {
      std::ofstream ofs(file_name,std::ios::app);
      ofs<<"# Coordinate "<<step<<std::endl;
      for(size_t i=0,i_size=cas.size();i<i_size;++i)
      {
         ofs<<i<<" "<<cas.at(i)<<std::endl;
      }
      ofs<<std::endl;
      ofs<<"# Bond "<<step<<std::endl;
      for(size_t i=1,i_size=cas.size();i<i_size;++i)
      {
         ofs<<i-1<<" "<<i<<" "<<step<<std::endl;
      }
      ofs<<std::endl<<std::endl;
      ofs.close();
   }
}
