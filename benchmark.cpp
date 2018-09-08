

#include "benchmark.h"

#include <fstream>

#include "lambdatwist/lambdatwist.p3p.h"
#include "kneip/kneip.h"
#include "utils/mlibtime.h"
#include <ke/ke.h>
#include <p3p_generator.h>
#include <ke/ke_utils.h>
#include <data.h>
#include <kneip/kneip_utils.h>
#include <mnumerics.h>
#include "utils/string_helpers.h"

/*********
 *
 *
 *
 * Standardize input output format,
 * The speed comparison is primarily between two methods that output rotation matrixes. Therefore use that as basis.
 * The ke uses a retarded format to maximize their speed...
 * convert to something resonable on output since it would always be followed by that anyways.
 *
 * Internally the methods should test the feasibility of the returned solutions, given the three points.
 * Ie dont return invalid solutions, nans, non rotations, transforms that move points behind cameras, or large reprojection errors
 * If they fail to do this it will always be done directly after.
 *
 * It would be very natural to use inheritance to create a nice modular structure for the benchmark. Doing so would hide the performance however,
 *  as vtable lookups take quite some time compared to the 300ns or so that the routines take to run... yes I tried it...
 *
 *
 * **/
#include <sstream>
#include <solver.h>



using std::cout;using std::endl;
using namespace cvl;

namespace mlib{
namespace klas{









template<class T,uint Rows,uint Cols>
std::string displayMatrix(cvl::Matrix<T,Rows,Cols> M,bool displayrowcol){

    std::vector<std::string> headers,rownames;
    if(displayrowcol){
        for(uint i=0;i<Cols;++i)
            headers.push_back(mlib::toStr(i));
        for(uint i=0;i<Rows;++i)
            rownames.push_back(mlib::toStr(i)+": ");
    }
    std::vector<std::vector<T>> rows;
    for(uint r=0;r<Rows;++r){
        std::vector<T> row;
        for(uint c=0;c<Cols;++c)
            row.push_back(M(r,c));
        rows.push_back(row);
    }

    return mlib::displayTable(headers,rows,rownames);
}





template<class T>class KeSolver{
public:
    int solve(Data<T>& data, cvl::Vector<Matrix<T,3,3>,4>& Rs,
              cvl::Vector<Vector<T,3>,4>& Ts){
        return kes::ke_p3p_fair(data,Rs,Ts);
    }
    std::string get_name(){return "Ke";}
};
template<class T>class KneipSolver{
public:
    int solve(Data<T>& data, cvl::Vector<Matrix<T,3,3>,4>& Rs,
              cvl::Vector<Vector<T,3>,4>& Ts){
        return kneip::kneip_p3p_fair(data,Rs,Ts);
    }
    std::string get_name(){return "Kneip";}
};

template<class T>class LambdaSolver{
public:
    int solve(Data<T>& data, cvl::Vector<Matrix<T,3,3>,4>& Rs,
              cvl::Vector<Vector<T,3>,4>& Ts){

        Vector3<Vector3<T>> yss=data.xr;
        Vector3<Vector3<T>> xss=data.x0;
        return p3p_lambdatwist(yss[0],yss[1],yss[2],xss[0],xss[1],xss[2],Rs,Ts);

    }
    std::string get_name(){return "Lambda";}
};













// verification of answers can be made slowly, so...

template<class T, class Solver> P3PResult<T> compute_accuracy(Solver S, std::vector<Data<T>> datas, T error_limit=1e-6){
    P3PResult<T> res(S.get_name());res.errors.reserve(datas.size());

    for(Data<T>& data:datas){
        cvl::Vector<Matrix<T,3,3>,4> Rs;
        cvl::Vector<Vector<T,3>,4> Ts;

        int valid = S.solve(data,Rs,Ts);
        res.valid+=valid; // valid according to solver.


        int duplicates=0;
        int sols=data.good_solutions(Rs,Ts,valid,duplicates); // valid according to data(harsher)
        if(sols==0)
            res.ground_truth_in_set++;
        if(valid>sols)
            res.incorrect_valid+=valid-sols;

        res.solutions+=sols;
        res.duplicates+=duplicates;

        T error=data.min_error(Rs,Ts,valid);
        res.errors.push_back(error);
        if(error<error_limit)
            res.ground_truth_in_set++;

    }
    return res;
}





template<class T>  void test(const std::vector<Data<T>>& datas){

    cout<<"Beginning Test: "<<endl;
    T error_limit=1e-6;


    P3PResult<T> lambda = compute_accuracy<T,LambdaSolver<T>>(LambdaSolver<T>(),datas,error_limit);
    P3PResult<T> ke = compute_accuracy<T,KeSolver<T>>(KeSolver<T>(),datas,error_limit);
    P3PResult<T> kneip= compute_accuracy<T,KneipSolver<T>>(KneipSolver<T>(),datas,error_limit);



    // make the result table
    std::vector<P3PResult<T>> res={lambda,ke,kneip};
    std::vector<std::string> headers;
    std::vector<std::string> columns;
    std::vector<std::vector<std::string>> rows;


    for(auto r:res) headers.push_back(r.name);

    {
        // number of times ground truth was found
        std::vector<int> gt;
        for(auto r:res) gt.push_back(r.ground_truth_in_set);
        rows.push_back(toStrVec(gt));
        columns.push_back("ground truth found");
    }
    {
        // ratio for ground truth found
        std::vector<double> gtratio;for(auto r:res) gtratio.push_back((double)r.ground_truth_in_set/datas.size());
        rows.push_back(toStrVec(gtratio));
        columns.push_back("ground truth found ratio");
    }
    {
        // number of no solution at all found
        std::vector<int> nosol;
        for(auto r:res) nosol.push_back(r.no_solution);
        rows.push_back(toStrVec(nosol));
        columns.push_back("no solution found");
    }
    {
        // number of valid solutions according to solver
        std::vector<int> solutions;
        for(auto r:res) solutions.push_back(r.valid);
        rows.push_back(toStrVec(solutions));
        columns.push_back("solutions");
    }
    {
        // duplicates
        std::vector<int> duplicates;
        for(auto r:res) duplicates.push_back(r.duplicates);
        rows.push_back(toStrVec(duplicates));
        columns.push_back("duplicates");
        //
    }
    {
        // duplicates
        std::vector<int> duplicates;
        for(auto r:res) duplicates.push_back(r.solutions - r.duplicates);
        rows.push_back(toStrVec(duplicates));
        columns.push_back("solutions - duplicates");
        //
    }
    {
        // incorrect valid solutions
        std::vector<int> solutions;
        for(auto r:res) solutions.push_back(r.incorrect_valid);
        rows.push_back(toStrVec(solutions));
        columns.push_back("incorrect valid solutions");
    }
    {
        // correct valid unique solutions
        std::vector<int> solutions;
        for(auto r:res) solutions.push_back(r.solutions - r.duplicates -r.incorrect_valid);
        rows.push_back(toStrVec(solutions));
        columns.push_back("correct valid unique");
    }


    cout<<displayTable(headers,rows,columns,"Table: ")<<endl;
    {
        std::vector<std::vector<T>> errors;
        for(auto r:res) errors.push_back(r.errors);
        std::string str=display(errors,headers);
    }

/*
    {
        // clear the screen
        cout << "\033[2J\033[1;1H"; cout.flush();


        std::vector<std::string> headers={"kneip","lambda","ke"};

        std::vector<std::vector<std::string>> rows;

        rows.push_back(toStrVec<int>({kneip.ground_truth_in_set, lambda.ground_truth_in_set, ke.ground_truth_in_set}));
        rows.push_back(toStrVec<double>({ (double)kneip.ground_truth_in_set/datas.size(), (double)lambda.ground_truth_in_set/datas.size(),(double)ke.ground_truth_in_set/datas.size()}));

        rows.push_back(toStrVec<int>({ kneip.valid, lambda.valid,ke.valid}));
        rows.push_back(toStrVec<int>({ kneip.get_incorrect_valid(),
                                       lambda.get_incorrect_valid(),
                                       ke.get_incorrect_valid()}));

        rows.push_back(toStrVec<int>({ kneip.solutions, lambda.solutions,ke.solutions}));
        rows.push_back(toStrVec<int>({ kneip.duplicates, lambda.duplicates,ke.duplicates}));

        rows.push_back(toStrVec<int>({kneip.solutions -kneip.duplicates, lambda.solutions -lambda.duplicates,ke.solutions -ke.duplicates}));

        cout<<displayTable(headers,rows,{"good: ","ratio: ","alg sols"," incorrect_sols: ", " verified solutions: ","dups:","unique correct solutions"},"Table: ")<<endl;

        std::vector<std::vector<T>> tmp={kneip.errors,lambda.errors,ke.errors};

        std::string str=display(tmp,{"kneip","lambda","ke"});

        cout<<str<<endl;
    }

    if(datas.size()<=10000000 && false){
        {
            {

                cout<<kneip.errors.size()<<endl;
                std::ofstream out("numerics.csv");
                //out<<"# row order is kneip,lambda, Ke\n";

                //  out<<mlib::getCSV<double>(kneip.errors)<<"\n";
                //  out<<mlib::getCSV<double>(lambda.errors)<<"\n";
                //  out<<mlib::getCSV<double>(ke.errors)<<"\n";
                out<<endl;
            }
        }
    }
    cout<<"Test - done! "<<endl;
    */
}






template<class T> void versus(const std::vector<Data<T>>& datas, int repeat_versus=5, bool inner_timers=false){

    cout<<"Beginning Versus"<<endl;


    Timer kneip("Kneip: ",datas.size());

    Timer lt("Lambda: ",datas.size());
    Timer kes("ke ",datas.size());

    Timer tot_lambda("Lambda Total");
    Timer tot_kes("ke Total");
    Timer tot_kneip("Kneip Total");



    // dont opt away
    int sols=0;
    // ke

    bool dokneip=true;
    for(int i=0;i<repeat_versus;++i){

        //kneip
        if(dokneip){

            tot_kneip.tic();
            for(int i=0;i<datas.size();++i){
                //if(inner_timers){if(i==100){                kneip.clear();            }kneip.tic();}

                cvl::Vector<Matrix<T,3,3>,4> Rs;
                cvl::Vector<Vector<T,3>,4> Ts;
                sols+= kneip::kneip_p3p_fair(datas[i],Rs, Ts);

                //    if(inner_timers)kneip.toc();
            }

            tot_kneip.toc();
        }


        {
            tot_kes.tic();

            for(int i=0;i<datas.size();++i){
                //    if(inner_timers)if(i==100){            kes.clear();        }        kes.tic();
                cvl::Vector<Matrix<T,3,3>,4> Rs;
                cvl::Vector<Vector<T,3>,4> Ts;
                sols+=kes::ke_p3p_fair(datas[i],Rs,Ts);
                //      if(inner_timers)kes.toc();

            }
            tot_kes.toc();
        }
        {
            tot_lambda.tic();
            for(int i=0;i<datas.size();++i){
                //    if(inner_timers) {if(i==100){                lt.clear();            }            lt.tic();}
                cvl::Vector<Matrix<T,3,3>,4> Rs;
                cvl::Vector<Vector<T,3>,4> Ts;
                Data<T> data=datas[i];
                Vector3<Vector<T,3>> yss=data.xr;
                Vector3<Vector<T,3>> xss=data.x0;
                sols+=p3p_lambdatwist(yss[0],yss[1],yss[2],xss[0],xss[1],xss[2],Rs,Ts);
                //  if(inner_timers)lt.toc();
            }
            tot_lambda.toc();
        }

        std::vector<mlib::Timer> tots={tot_lambda,tot_kes,tot_kneip};

        printtimers();
        std::vector<mlib::Timer> ts={lt,kneip,kes};
        cout<<tots<<endl<<"wrote the ts"<<endl;
    }
    if(false){

        std::ofstream out("timing.csv");

        {
            auto ts=kneip.getTimes();
            for(Time time:ts){

                out<<(int)(time.ns)<<", ";
            }
            out<<"\n";
        }
        {
            auto ts=lt.getTimes();
            for(Time time:ts){

                out<<(int)(time.ns)<<", ";
            }
            out<<"\n";
        }
        {
            auto ts=kes.getTimes();
            for(Time time:ts){

                out<<(int)(time.ns)<<", ";
            }
            out<<"\n";
        }
        out<<std::endl;
        out.close();
        //cout<<"means: "<<mean(kn)<<" "<<mean(ls)<<endl;

    }

    //cout<<mlib::display(samples,false)<<endl;
    cout<<"Versus Done !"<<endl;
}



template<class T> void test_type(){
    // generate data
    cout<<"generating";cout.flush();
    Generator gennie;

    int N=std::pow(10,7);

    std::vector<Data<T>> datas;datas.reserve(N);
    for(int i=0;i<N;++i){
        datas.push_back(gennie.next<T>());
        //if((i%1000)==999)        cout<<"i: "<<i<<endl;
    }
    cout<<" done "<<datas.size()<<endl;

    test(datas);
    // printtimers();

    versus(datas);

    //testReal();
}


template<class T> void profile_lambda(){

    Generator gennie;

    int sols=0;
    cvl::Vector<Matrix<T,3,3>,4> Rs;
    cvl::Vector<Vector<T,3>,4> Ts;
    Data<T> data=gennie.next<T>();

    for(int i=0;i<100000;++i){

        Vector3<Vector3<T>> yss=data.xr;
        Vector3<Vector3<T>> xss=data.x0;
        sols+=p3p_lambdatwist(yss[0],yss[1],yss[2],xss[0],xss[1],xss[2],Rs,Ts);
    }
    cout<<sols<<endl;
}




void testAll(){


    test_type<double>();


}
}// end namespace klas
}// end namespace mlib
