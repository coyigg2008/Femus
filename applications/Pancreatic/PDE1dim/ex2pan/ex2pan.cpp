//Solve - \Delta u = 1
#include "FemusInit.hpp"
#include "MultiLevelSolution.hpp"
#include "MultiLevelProblem.hpp"
#include "NumericVector.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "TransientSystem.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "adept.h"
// #include "CurrentElem.hpp"
// #include "ElemType_template.hpp"
#include "pancreatic_parameter.hpp"

using namespace std;
using namespace femus;

double InitialValueC(const std::vector < double >& x) {
  return 0.4;
}
double InitialValueT4(const std::vector < double >& x) {
  return 1.0e-3;
}
double InitialValueT8(const std::vector < double >& x) {
  return 1.05e-3;
}
double InitialValueTr(const std::vector < double >& x) {
  return 9.0e-4;
}
double InitialValueD(const std::vector < double >& x) {
  return 5.0e-5;
}
double InitialValueR(const std::vector < double >& x) {
  return 0.01;
}

bool SetBoundaryCondition(const std::vector < double >& x, const char name[], double& value, const int face_name, const double time) {

  bool dirichlet = false;
  value = 0;
  
//   const double tolerance = 1.e-5;
//   
//   if (face_name == 1) {
//       dirichlet = true;
//         value = 0.;
//   }
  
  return dirichlet;
 }

void AssembleBoussinesqAppoximation_AD(MultiLevelProblem& ml_prob);

int main(int argc, char** args) {

  // init Petsc-MPI communicator
  FemusInit mpinit(argc, args, MPI_COMM_WORLD);
  MultiLevelMesh ml_mesh;
  
  std::string fe_quad_rule("seventh");
  ml_mesh.GenerateCoarseBoxMesh(2,2,0,0.,4.,0.,4.,0.,0.,QUAD9,fe_quad_rule.c_str());
  // ml_mesh.GenerateCoarseBoxMesh(2,2,0,0.,1.,0.,1.,0.,0.,QUAD9,fe_quad_rule.c_str());
  //ml_mesh.GenerateCoarseBoxMesh(2,2,2,0.,1.,0.,1.,0.,1.,HEX27,fe_quad_rule.c_str());
  unsigned dim = ml_mesh.GetDimension();
  
  unsigned numberOfUniformLevels = 3;
  unsigned numberOfSelectiveLevels = 0;
  ml_mesh.RefineMesh(numberOfUniformLevels , numberOfUniformLevels + numberOfSelectiveLevels, NULL);
// ml_mesh.EraseCoarseLevels(3);
//  ml_mesh.EraseCoarseLevels(numberOfUniformLevels + numberOfSelectiveLevels - 1);
  
  //print mesh info
  ml_mesh.PrintInfo();
  MultiLevelSolution mlSol(&ml_mesh);
  
  // add variables to mlSol
  mlSol.AddSolution("C", LAGRANGE, SECOND,2);
  mlSol.AddSolution("D", LAGRANGE, SECOND,2);
  mlSol.AddSolution("T4", LAGRANGE, SECOND,2);
  mlSol.AddSolution("T8", LAGRANGE, SECOND,2);
  mlSol.AddSolution("Tr", LAGRANGE, SECOND,2);
  mlSol.AddSolution("I2", LAGRANGE, SECOND,2);
  mlSol.AddSolution("I12", LAGRANGE, SECOND,2); 
  mlSol.AddSolution("EC", LAGRANGE, SECOND,2); 
  mlSol.AddSolution("ED", LAGRANGE, SECOND,2);
  mlSol.AddSolution("Talpha", LAGRANGE, SECOND,2);
  mlSol.AddSolution("m1", LAGRANGE, SECOND,2);
  mlSol.AddSolution("m2", LAGRANGE, SECOND,2);
  mlSol.AddSolution("S", LAGRANGE, SECOND,2);
  mlSol.AddSolution("u", LAGRANGE, SECOND,2);
  mlSol.AddSolution("R", LAGRANGE, SECOND,2);
  
  mlSol.Initialize("All");
  mlSol.Initialize("C", InitialValueC); 
  mlSol.Initialize("T4",InitialValueT4);
  mlSol.Initialize("T8",InitialValueT8);
  mlSol.Initialize("Tr",InitialValueTr);
  mlSol.Initialize("D", InitialValueD);
  
  // attach the boundary condition function and generate boundary data
  mlSol.AttachSetBoundaryConditionFunction(SetBoundaryCondition);
  mlSol.GenerateBdc("All");
//   mlSol.GenerateBdc("C");
//   mlSol.GenerateBdc("D");
  
  // define the multilevel problem attach the mlSol object to it
  MultiLevelProblem mlProb(&mlSol);
  
    // add system Poisson in mlProb as a Linear Implicit System
  TransientNonlinearImplicitSystem& system = mlProb.add_system < TransientNonlinearImplicitSystem > ("NS");
  
  system.AddSolutionToSystemPDE("C");
  system.AddSolutionToSystemPDE("D");
  system.AddSolutionToSystemPDE("T4");
  system.AddSolutionToSystemPDE("T8");
  system.AddSolutionToSystemPDE("Tr");
  system.AddSolutionToSystemPDE("I2");
  system.AddSolutionToSystemPDE("I12");
  system.AddSolutionToSystemPDE("EC");
  system.AddSolutionToSystemPDE("ED");
  system.AddSolutionToSystemPDE("Talpha");
  system.AddSolutionToSystemPDE("m1");
  system.AddSolutionToSystemPDE("m2");
  system.AddSolutionToSystemPDE("S");
  system.AddSolutionToSystemPDE("u");
  system.AddSolutionToSystemPDE("R");
  
  system.SetLinearEquationSolverType(FEMuS_DEFAULT);
  // attach the assembling function to system
  system.SetAssembleFunction(AssembleBoussinesqAppoximation_AD);
  
  system.SetMaxNumberOfNonLinearIterations(20);
  system.SetNonLinearConvergenceTolerance(1.e-6);
  
  system.SetMaxNumberOfLinearIterations(1);
  system.SetAbsoluteLinearConvergenceTolerance(1.e-15);
  
  system.SetMgType(F_CYCLE);

  system.SetNumberPreSmoothingStep(1);
  system.SetNumberPostSmoothingStep(1);
  //initialize
  system.init();
  
  system.SetSolverFineGrids(RICHARDSON);
  system.SetPreconditionerFineGrids(ILU_PRECOND);
//  system.SetRichardsonScaleFactor(.6);
  system.SetTolerances(1.e-8, 1.e-10, 1.e+50, 100, 100);
  
  system.ClearVariablesToBeSolved();
  system.AddVariableToBeSolved("All");
//  system.SetElementBlockNumber("All");
//   system.SetNumberOfSchurVariables(1);
  
//   std::vector< double > x(3);
//   x[0] = 1.5; //the marker is in element 117 (proc 1)
//   x[1] = 0.0;
//   x[2] = 0.0;
  
  // print solutions
  std::vector < std::string > variablesToBePrinted;
  variablesToBePrinted.push_back("All");

  VTKWriter vtkIO(&mlSol);
  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, 0);
  
  double dt = 0.0005;
  system.SetIntervalTime(dt);
  unsigned n_timesteps = 80000;
  
//   Marker marker(x, 0., VOLUME, mlSol.GetLevel(numberOfUniformLevels - 1), dim , true); % change the 2 to dim, becuase I believe it is one dimension
//   unsigned elem = marker.GetMarkerElement();
//   std::vector<double> xi;
//   marker.GetMarkerLocalCoordinates(xi);
//   
//   
//   char out_file1[100]="";
//   strcpy(out_file1,"C.dat");
//   ofstream outfile1(out_file1,ios::out|ios::trunc|ios::binary);
// 
//   char out_file2[100]="";
//   strcpy(out_file2,"D.dat");
//   ofstream outfile2(out_file2,ios::out|ios::trunc|ios::binary);
//   
//   char out_file3[100]="";
//   strcpy(out_file3,"T4.dat");
//   ofstream outfile3(out_file3,ios::out|ios::trunc|ios::binary);
//   
//   char out_file4[100]="";
//   strcpy(out_file4,"T8.dat");
//   ofstream outfile4(out_file4,ios::out|ios::trunc|ios::binary);
//   
//   char out_file5[100]="";
//   strcpy(out_file5,"Tr.dat");
//   ofstream outfile5(out_file5,ios::out|ios::trunc|ios::binary);
//   
//   char out_file6[100]="";
//   strcpy(out_file6,"I2.dat");
//   ofstream outfile6(out_file6,ios::out|ios::trunc|ios::binary);
//   
//   char out_file7[100]="";
//   strcpy(out_file7,"I12.dat");
//   ofstream outfile7(out_file7,ios::out|ios::trunc|ios::binary);
//   
//   char out_file8[100]="";
//   strcpy(out_file8,"EC.dat");
//   ofstream outfile8(out_file8,ios::out|ios::trunc|ios::binary);
//   
//   char out_file9[100]="";
//   strcpy(out_file9,"ED.dat");
//   ofstream outfile9(out_file9,ios::out|ios::trunc|ios::binary);
//   
//   char out_file10[100]="";
//   strcpy(out_file10,"Talpha.dat");
//   ofstream outfile10(out_file10,ios::out|ios::trunc|ios::binary);
//   
//   char out_file11[100]="";
//   strcpy(out_file11,"m1.dat");
//   ofstream outfile11(out_file11,ios::out|ios::trunc|ios::binary);
//   
//   char out_file12[100]="";
//   strcpy(out_file12,"m2.dat");
//   ofstream outfile12(out_file12,ios::out|ios::trunc|ios::binary);
//   
//   char out_file13[100]="";
//   strcpy(out_file13,"S.dat");
//   ofstream outfile13(out_file13,ios::out|ios::trunc|ios::binary);
//   
//   char out_file14[100]="";
//   strcpy(out_file14,"u.dat");
//   ofstream outfile14(out_file14,ios::out|ios::trunc|ios::binary);
//   
//   char out_file15[100]="";
//   strcpy(out_file15,"R.dat");
//   ofstream outfile15(out_file15,ios::out|ios::trunc|ios::binary);
//   
//   vector <double> sol_Fpara(8);
//   vector <double> sol_Spara(7);
//   std::pair < vector <double>, vector <double> > out_value;
  
  for(unsigned time_step = 0; time_step < n_timesteps; time_step++) {
 
     if(time_step > 0) system.SetMgType(V_CYCLE);
 
     system.MGsolve();
     system.CopySolutionToOldSolution();
//      out_value = GetVaribleValues(mlProb, elem, xi);
//      sol_Fpara = out_value.first;
//      sol_Spara = out_value.second;
//      
//      outfile1 << (time_step + 1) * dt <<"  "<< sol_Fpara[0] << std::endl;
//      outfile2 << (time_step + 1) * dt <<"  "<< sol_Fpara[1] << std::endl;
//      outfile3 << (time_step + 1) * dt <<"  "<< sol_Fpara[2] << std::endl;
//      outfile4 << (time_step + 1) * dt <<"  "<< sol_Fpara[3] << std::endl;
//      outfile5 << (time_step + 1) * dt <<"  "<< sol_Fpara[4] << std::endl;
//      outfile6 << (time_step + 1) * dt <<"  "<< sol_Fpara[5] << std::endl;
//      outfile7 << (time_step + 1) * dt <<"  "<< sol_Fpara[6] << std::endl;
//      outfile8 << (time_step + 1) * dt <<"  "<< sol_Fpara[7] << std::endl;
//      
//      outfile9  << (time_step + 1) * dt <<"  "<< sol_Spara[0] << std::endl;
//      outfile10 << (time_step + 1) * dt <<"  "<< sol_Spara[1] << std::endl;
//      outfile11 << (time_step + 1) * dt <<"  "<< sol_Spara[2] << std::endl;
//      outfile12 << (time_step + 1) * dt <<"  "<< sol_Spara[3] << std::endl;
//      outfile13 << (time_step + 1) * dt <<"  "<< sol_Spara[4] << std::endl;
//      outfile14 << (time_step + 1) * dt <<"  "<< sol_Spara[5] << std::endl;
//      outfile15 << (time_step + 1) * dt <<"  "<< sol_Spara[6] << std::endl;
     if ((time_step + 1) % 1000 ==0)  vtkIO.Write(DEFAULT_OUTPUTDIR, "biquadratic", variablesToBePrinted, time_step + 1);
   }
//   outfile1.close();
//   outfile2.close();
//   outfile3.close();
//   outfile4.close();
//   outfile5.close();
//   outfile6.close();
//   outfile7.close();
//   outfile8.close();
//   outfile9.close();
//   outfile10.close();
//   outfile11.close();
//   outfile12.close();
//   outfile13.close();
//   outfile14.close();
//   outfile15.close();
  ml_mesh.PrintInfo();
  
  return 0;
}


void AssembleBoussinesqAppoximation_AD(MultiLevelProblem& ml_prob) {
  //  ml_prob is the global object from/to where get/set all the data
  //  level is the level of the PDE system to be assembled
  //  levelMax is the Maximum level of the MultiLevelProblem
  //  assembleMatrix is a flag that tells if only the residual or also the matrix should be assembled

  //  extract pointers to the several objects that we are going to use
  TransientNonlinearImplicitSystem* mlPdeSys   = &ml_prob.get_system<TransientNonlinearImplicitSystem> ("NS");   // pointer to the linear implicit system named "Poisson"
  const unsigned level = mlPdeSys->GetLevelToAssemble();

  Mesh*           msh         = ml_prob._ml_msh->GetLevel(level);    // pointer to the mesh (level) object
  elem*           el          = msh->el;  // pointer to the elem object in msh (level)

  MultiLevelSolution*   mlSol         = ml_prob._ml_sol;  // pointer to the multilevel solution object
  Solution*   sol         = ml_prob._ml_sol->GetSolutionLevel(level);    // pointer to the solution (level) object


  LinearEquationSolver* pdeSys        = mlPdeSys->_LinSolver[level];  // pointer to the equation (level) object

  bool assembleMatrix = mlPdeSys->GetAssembleMatrix();
  // call the adept stack object
  adept::Stack& s = FemusInit::_adeptStack;
  if(assembleMatrix) s.continue_recording();
  else s.pause_recording();

  SparseMatrix*   KK          = pdeSys->_KK;  // pointer to the global stifness matrix object in pdeSys (level)
  NumericVector*  RES         = pdeSys->_RES; // pointer to the global residual vector object in pdeSys (level)

  const unsigned  dim = msh->GetDimension(); // get the domain dimension of the problem
  unsigned dim2 = (3 * (dim - 1) + !(dim - 1));        // dim2 is the number of second order partial derivatives (1,3,6 depending on the dimension)
  unsigned    iproc = msh->processor_id(); // get the process_id (for parallel computation)

  // reserve memory for the local standar vectors
  const unsigned maxSize = static_cast< unsigned >(ceil(pow(3, dim)));          // conservative: based on line3, quad9, hex27

  unsigned solCIndex, solDIndex, solT4Index, solT8Index, solTrIndex, 
           solI2Index, solI12Index, solECIndex, solEDIndex, solTalphaIndex, 
           solm1Index, solm2Index, solSIndex, soluIndex, solRIndex;     
  
  solCIndex         = mlSol->GetIndex("C");    // get the position of "C" in the ml_sol object
  solDIndex         = mlSol->GetIndex("D");
  solT4Index        = mlSol->GetIndex("T4");
  solT8Index        = mlSol->GetIndex("T8");
  solTrIndex        = mlSol->GetIndex("Tr");
  solI2Index        = mlSol->GetIndex("I2");
  solI12Index       = mlSol->GetIndex("I12");
  solECIndex        = mlSol->GetIndex("EC");
  solEDIndex        = mlSol->GetIndex("ED");
  solTalphaIndex    = mlSol->GetIndex("Talpha");
  solm1Index        = mlSol->GetIndex("m1");
  solm2Index        = mlSol->GetIndex("m2");
  solSIndex         = mlSol->GetIndex("S");
  soluIndex         = mlSol->GetIndex("u");
  solRIndex         = mlSol->GetIndex("R");
  
  unsigned solCType      = mlSol->GetSolutionType(solCIndex);    // get the finite element type for "T"
  unsigned solDType      = mlSol->GetSolutionType(solDIndex);
  unsigned solT4Type     = mlSol->GetSolutionType(solT4Index);
  unsigned solT8Type     = mlSol->GetSolutionType(solT8Index);
  unsigned solTrType     = mlSol->GetSolutionType(solTrIndex);
  unsigned solI2Type     = mlSol->GetSolutionType(solI2Index);
  unsigned solI12Type    = mlSol->GetSolutionType(solI12Index);
  unsigned solECType     = mlSol->GetSolutionType(solECIndex);
  unsigned solEDType     = mlSol->GetSolutionType(solEDIndex);
  unsigned solTalphaType = mlSol->GetSolutionType(solTalphaIndex);
  unsigned solm1Type     = mlSol->GetSolutionType(solm1Index);
  unsigned solm2Type     = mlSol->GetSolutionType(solm2Index);
  unsigned solSType      = mlSol->GetSolutionType(solSIndex);
  unsigned soluType      = mlSol->GetSolutionType(soluIndex);
  unsigned solRType      = mlSol->GetSolutionType(solRIndex);
  
  unsigned solCPdeIndex,      solDPdeIndex,    solT4PdeIndex, solT8PdeIndex, solTrPdeIndex, 
           solI2PdeIndex,     solI12PdeIndex,  solECPdeIndex, solEDPdeIndex, solTalphaPdeIndex, 
           solm1PdeIndex,     solm2PdeIndex,   solSPdeIndex,  soluPdeIndex,  solRPdeIndex; 
  
 solCPdeIndex       = mlPdeSys->GetSolPdeIndex("C");    // get the position of "T" in the pdeSys object
 solDPdeIndex       = mlPdeSys->GetSolPdeIndex("D");
 solT4PdeIndex      = mlPdeSys->GetSolPdeIndex("T4");
 solT8PdeIndex      = mlPdeSys->GetSolPdeIndex("T8");
 solTrPdeIndex      = mlPdeSys->GetSolPdeIndex("Tr");
 solI2PdeIndex      = mlPdeSys->GetSolPdeIndex("I2");
 solI12PdeIndex     = mlPdeSys->GetSolPdeIndex("I12");
 solECPdeIndex      = mlPdeSys->GetSolPdeIndex("EC");
 solEDPdeIndex      = mlPdeSys->GetSolPdeIndex("ED");
 solTalphaPdeIndex  = mlPdeSys->GetSolPdeIndex("Talpha");
 solm1PdeIndex      = mlPdeSys->GetSolPdeIndex("m1");
 solm2PdeIndex      = mlPdeSys->GetSolPdeIndex("m2");
 solSPdeIndex       = mlPdeSys->GetSolPdeIndex("S");
 soluPdeIndex       = mlPdeSys->GetSolPdeIndex("u");
 solRPdeIndex       = mlPdeSys->GetSolPdeIndex("R");

  vector < adept::adouble >  solC; // local solution
  vector < adept::adouble >  solD;
  vector < adept::adouble >  solT4;
  vector < adept::adouble >  solT8;
  vector < adept::adouble >  solTr;
  vector < adept::adouble >  solI2;
  vector < adept::adouble >  solI12;
  vector < adept::adouble >  solEC;
  vector < adept::adouble >  solED;
  vector < adept::adouble >  solTalpha;
  vector < adept::adouble >  solm1;
  vector < adept::adouble >  solm2;
  vector < adept::adouble >  solS;
  vector < adept::adouble >  solu;
  vector < adept::adouble >  solR;
  
  vector < double >  solCold; // local old solution
  vector < double >  solDold;
  vector < double >  solT4old;
  vector < double >  solT8old;
  vector < double >  solTrold;
  vector < double >  solI2old;
  vector < double >  solI12old;
  vector < double >  solECold;
  vector < double >  solEDold;
  vector < double >  solTalphaold;
  vector < double >  solm1old;
  vector < double >  solm2old;
  vector < double >  solSold;
  vector < double >  soluold;
  vector < double >  solRold;
  
  vector < adept::adouble >  aResC; // local redidual vector
  vector < adept::adouble >  aResD;
  vector < adept::adouble >  aResT4;
  vector < adept::adouble >  aResT8;
  vector < adept::adouble >  aResTr;
  vector < adept::adouble >  aResI2;
  vector < adept::adouble >  aResI12;
  vector < adept::adouble >  aResEC;
  vector < adept::adouble >  aResED;
  vector < adept::adouble >  aResTalpha;
  vector < adept::adouble >  aResm1;
  vector < adept::adouble >  aResm2;
  vector < adept::adouble >  aResS;
  vector < adept::adouble >  aResu;
  vector < adept::adouble >  aResR;

  vector < vector < double > > coordX(dim);    // local coordinates
  unsigned coordXType = 2; // get the finite element type for "x", it is always 2 (LAGRANGE QUADRATIC)

  solC.reserve(maxSize); //allocate the meomery 
  solD.reserve(maxSize);
  solT4.reserve(maxSize);
  solT8.reserve(maxSize);
  solTr.reserve(maxSize);
  solI2.reserve(maxSize);
  solI12.reserve(maxSize);
  solEC.reserve(maxSize);
  solED.reserve(maxSize);
  solTalpha.reserve(maxSize);
  solm1.reserve(maxSize);
  solm2.reserve(maxSize);
  solS.reserve(maxSize);
  solu.reserve(maxSize);
  solR.reserve(maxSize);
  
  solCold.reserve(maxSize);
  solDold.reserve(maxSize);
  solT4old.reserve(maxSize);
  solT8old.reserve(maxSize);
  solTrold.reserve(maxSize);
  solI2old.reserve(maxSize);
  solI12old.reserve(maxSize);
  solECold.reserve(maxSize);
  solEDold.reserve(maxSize);
  solTalphaold.reserve(maxSize);
  solm1old.reserve(maxSize);
  solm2old.reserve(maxSize);
  solSold.reserve(maxSize);
  soluold.reserve(maxSize);
  solRold.reserve(maxSize);
  
  aResC.reserve(maxSize);
  aResD.reserve(maxSize);
  aResT4.reserve(maxSize);
  aResT8.reserve(maxSize);
  aResTr.reserve(maxSize);
  aResI2.reserve(maxSize);
  aResI12.reserve(maxSize);
  aResEC.reserve(maxSize);
  aResED.reserve(maxSize);
  aResTalpha.reserve(maxSize);
  aResm1.reserve(maxSize);
  aResm2.reserve(maxSize);
  aResS.reserve(maxSize);
  aResu.reserve(maxSize);
  aResR.reserve(maxSize);

  vector <double> phiC;  // local test function
  vector <double> phiC_x; // local test function first order partial derivatives
  vector <double> phiC_xx; // local test function second order partial derivatives
  
  vector <double> phiD;  // local test function
  vector <double> phiD_x; // local test function first order partial derivatives
  vector <double> phiD_xx; // local test function second order partial derivatives
  
  vector <double> phiT4;  // local test function
  vector <double> phiT4_x; // local test function first order partial derivatives
  vector <double> phiT4_xx; // local test function second order partial derivatives
  
  vector <double> phiT8;  // local test function
  vector <double> phiT8_x; // local test function first order partial derivatives
  vector <double> phiT8_xx; // local test function second order partial derivatives
  
  vector <double> phiTr;  // local test function
  vector <double> phiTr_x; // local test function first order partial derivatives
  vector <double> phiTr_xx; // local test function second order partial derivatives
  
  vector <double> phiI2;  // local test function
  vector <double> phiI2_x; // local test function first order partial derivatives
  vector <double> phiI2_xx; // local test function second order partial derivatives
  
  vector <double> phiI12;  // local test function
  vector <double> phiI12_x; // local test function first order partial derivatives
  vector <double> phiI12_xx; // local test function second order partial derivatives
  
  vector <double> phiEC;  // local test function
  vector <double> phiEC_x; // local test function first order partial derivatives
  vector <double> phiEC_xx; // local test function second order partial derivatives

  vector <double> phiED;  // local test function
  vector <double> phiED_x; // local test function first order partial derivatives
  vector <double> phiED_xx; // local test function second order partial derivatives
  
  vector <double> phiTalpha;  // local test function
  vector <double> phiTalpha_x; // local test function first order partial derivatives
  vector <double> phiTalpha_xx; // local test function second order partial derivatives

  vector <double> phim1;  // local test function
  vector <double> phim1_x; // local test function first order partial derivatives
  vector <double> phim1_xx; // local test function second order partial derivatives
  
  vector <double> phim2;  // local test function
  vector <double> phim2_x; // local test function first order partial derivatives
  vector <double> phim2_xx; // local test function second order partial derivatives
  
  vector <double> phiS;  // local test function
  vector <double> phiS_x; // local test function first order partial derivatives
  vector <double> phiS_xx; // local test function second order partial derivatives
  
  vector <double> phiu;  // local test function
  vector <double> phiu_x; // local test function first order partial derivatives
  vector <double> phiu_xx; // local test function second order partial derivatives
  
  vector <double> phiR;  // local test function
  vector <double> phiR_x; // local test function first order partial derivatives
  vector <double> phiR_xx; // local test function second order partial derivatives
  
  phiC.reserve(maxSize); // allocate the meomery
  phiC_x.reserve(maxSize * dim);
  phiC_xx.reserve(maxSize * dim2);
  
  phiD.reserve(maxSize); // allocate the meomery
  phiD_x.reserve(maxSize * dim);
  phiD_xx.reserve(maxSize * dim2);
  
  phiT4.reserve(maxSize); // allocate the meomery
  phiT4_x.reserve(maxSize * dim);
  phiT4_xx.reserve(maxSize * dim2);
  
  phiT8.reserve(maxSize); // allocate the meomery
  phiT8_x.reserve(maxSize * dim);
  phiT8_xx.reserve(maxSize * dim2);
  
  phiTr.reserve(maxSize); // allocate the meomery
  phiTr_x.reserve(maxSize * dim);
  phiTr_xx.reserve(maxSize * dim2);
  
  phiI2.reserve(maxSize); // allocate the meomery
  phiI2_x.reserve(maxSize * dim);
  phiI2_xx.reserve(maxSize * dim2);
  
  phiI12.reserve(maxSize); // allocate the meomery
  phiI12_x.reserve(maxSize * dim);
  phiI12_xx.reserve(maxSize * dim2);
  
  phiEC.reserve(maxSize); // allocate the meomery
  phiEC_x.reserve(maxSize * dim);
  phiEC_xx.reserve(maxSize * dim2);
  
  phiED.reserve(maxSize); // allocate the meomery
  phiED_x.reserve(maxSize * dim);
  phiED_xx.reserve(maxSize * dim2);
  
  phiTalpha.reserve(maxSize); // allocate the meomery
  phiTalpha_x.reserve(maxSize * dim);
  phiTalpha_xx.reserve(maxSize * dim2);

  phim1.reserve(maxSize); // allocate the meomery
  phim1_x.reserve(maxSize * dim);
  phim1_xx.reserve(maxSize * dim2);
  
  phim2.reserve(maxSize); // allocate the meomery
  phim2_x.reserve(maxSize * dim);
  phim2_xx.reserve(maxSize * dim2);
  
  phiS.reserve(maxSize); // allocate the meomery
  phiS_x.reserve(maxSize * dim);
  phiS_xx.reserve(maxSize * dim2);
  
  phiu.reserve(maxSize); // allocate the meomery
  phiu_x.reserve(maxSize * dim);
  phiu_xx.reserve(maxSize * dim2);
  
  phiR.reserve(maxSize); // allocate the meomery
  phiR_x.reserve(maxSize * dim);
  phiR_xx.reserve(maxSize * dim2);
  
  double weight; // gauss point weight

  vector< int > sysDof; // local to global pdeSys dofs
  sysDof.reserve(15 * maxSize);

  vector< double > Res; // local redidual vector
  Res.reserve(15 * maxSize);

  vector < double > Jac;
  Jac.reserve(15 * maxSize * 15 *maxSize);

  if(assembleMatrix) KK->zero(); // Set to zero all the entries of the Global Matrix

 
  // element loop: each process loops only on the elements that owns
  for(int iel = msh->_elementOffset[iproc]; iel < msh->_elementOffset[iproc + 1]; iel++) {
      
    // element geometry type
    short unsigned ielGeom = msh->GetElementType(iel);
    unsigned nDofsC = msh->GetElementDofNumber(iel, solCType);    // number of solution element dofs
    unsigned nDofsD = msh->GetElementDofNumber(iel, solDType);    // number of solution element dofs
    unsigned nDofsT4 = msh->GetElementDofNumber(iel, solT4Type);    // number of solution element dofs
    unsigned nDofsT8 = msh->GetElementDofNumber(iel, solT8Type);    // number of solution element dofs
    unsigned nDofsTr = msh->GetElementDofNumber(iel, solTrType);    // number of solution element dofs
    unsigned nDofsI2 = msh->GetElementDofNumber(iel, solI2Type);    // number of solution element dofs
    unsigned nDofsI12 = msh->GetElementDofNumber(iel, solI12Type);    // number of solution element dofs
    unsigned nDofsEC = msh->GetElementDofNumber(iel, solECType);    // number of solution element dofs
    unsigned nDofsED = msh->GetElementDofNumber(iel, solEDType);    // number of solution element dofs
    unsigned nDofsTalpha = msh->GetElementDofNumber(iel, solTalphaType);    // number of solution element dofs
    unsigned nDofsm1 = msh->GetElementDofNumber(iel, solm1Type);    // number of solution element dofs
    unsigned nDofsm2 = msh->GetElementDofNumber(iel, solm2Type);    // number of solution element dofs
    unsigned nDofsS = msh->GetElementDofNumber(iel, solSType);    // number of solution element dofs
    unsigned nDofsu = msh->GetElementDofNumber(iel, soluType);    // number of solution element dofs
    unsigned nDofsR = msh->GetElementDofNumber(iel, solRType);    // number of solution element dofs
    unsigned nDofsX = msh->GetElementDofNumber(iel, coordXType);    // number of coordinate element dofs

    unsigned nDofsTot = nDofsC + nDofsD + nDofsT4 + nDofsT8 + nDofsTr + nDofsI2 + nDofsI12 
                       + nDofsEC + nDofsED + nDofsTalpha + nDofsm1 + nDofsm2 + nDofsS + nDofsu + nDofsR;
    // resize local arrays
    sysDof.resize(nDofsTot);
    
    solC.resize(nDofsC);
    solD.resize(nDofsD);
    solT4.resize(nDofsT4);
    solT8.resize(nDofsT8);
    solTr.resize(nDofsTr);
    solI2.resize(nDofsI2);
    solI12.resize(nDofsI12);
    solEC.resize(nDofsEC);
    solED.resize(nDofsED);
    solTalpha.resize(nDofsTalpha);
    solm1.resize(nDofsm1);
    solm2.resize(nDofsm2);
    solS.resize(nDofsS);
    solu.resize(nDofsu);
    solR.resize(nDofsR);
    
    for(unsigned  k = 0; k < dim; k++) coordX[k].resize(nDofsX);
    solCold.resize(nDofsC);
    solDold.resize(nDofsD);
    solT4old.resize(nDofsT4);
    solT8old.resize(nDofsT8);
    solTrold.resize(nDofsTr);
    solI2old.resize(nDofsI2);
    solI12old.resize(nDofsI12);
    solECold.resize(nDofsEC);
    solEDold.resize(nDofsED);
    solTalphaold.resize(nDofsTalpha);
    solm1old.resize(nDofsm1);
    solm2old.resize(nDofsm2);
    solSold.resize(nDofsS);
    soluold.resize(nDofsu);
    solRold.resize(nDofsR);
    
    aResC.resize(nDofsC);
    aResD.resize(nDofsD);
    aResT4.resize(nDofsT4);
    aResT8.resize(nDofsT8);
    aResTr.resize(nDofsTr);
    aResI2.resize(nDofsI2);
    aResI12.resize(nDofsI12);
    aResEC.resize(nDofsEC);
    aResED.resize(nDofsED);
    aResTalpha.resize(nDofsTalpha);
    aResm1.resize(nDofsm1);
    aResm2.resize(nDofsm2);
    aResS.resize(nDofsS);
    aResu.resize(nDofsu);
    aResR.resize(nDofsR);
    
    std::fill(aResC.begin(),   aResC.end(), 0);    //set aRes to zero
    std::fill(aResD.begin(),   aResD.end(), 0);    //set aRes to zero
    std::fill(aResT4.begin(),  aResT4.end(), 0);    //set aRes to zero
    std::fill(aResT8.begin(),  aResT8.end(), 0);    //set aRes to zero
    std::fill(aResTr.begin(),  aResTr.end(), 0);    //set aRes to zero
    std::fill(aResI2.begin(),  aResI2.end(), 0);    //set aRes to zero
    std::fill(aResI12.begin(), aResI12.end(), 0);    //set aRes to zero
    std::fill(aResEC.begin(),  aResEC.end(), 0);    //set aRes to zero
    std::fill(aResED.begin(),  aResED.end(), 0);    //set aRes to zero
    std::fill(aResTalpha.begin(),  aResTalpha.end(), 0);    //set aRes to zero
    std::fill(aResm1.begin(),  aResm1.end(), 0);    //set aRes to zero
    std::fill(aResm2.begin(),  aResm2.end(), 0);    //set aRes to zero
    std::fill(aResS.begin(),   aResS.end(), 0);    //set aRes to zero
    std::fill(aResu.begin(),   aResu.end(), 0);    //set aRes to zero
    std::fill(aResR.begin(),   aResR.end(), 0);    //set aRes to zero
    
// local storage of global mapping and solution
    for(unsigned i = 0; i < nDofsC; i++) {
      unsigned solCDof = msh->GetSolutionDof(i, iel, solCType);    // global to global mapping between solution node and solution dof
      solC[i] = (*sol->_Sol[solCIndex])(solCDof);      // global extraction and local storage for the solution 
      solCold[i] = (*sol->_SolOld[solCIndex])(solCDof);      // global extraction and local storage for the solution
      sysDof[i] = pdeSys->GetSystemDof(solCIndex, solCPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }

    for(unsigned i = 0; i < nDofsD; i++) {
      unsigned solDDof = msh->GetSolutionDof(i, iel, solDType);    // global to global mapping between solution node and solution dof
      solD[i] = (*sol->_Sol[solDIndex])(solDDof);      // global extraction and local storage for the solution
      solDold[i] = (*sol->_SolOld[solDIndex])(solDDof);      // global extraction and local storage for the solution
      sysDof[i+nDofsC] = pdeSys->GetSystemDof(solDIndex, solDPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }

    for(unsigned i = 0; i < nDofsT4; i++) {
      unsigned solT4Dof = msh->GetSolutionDof(i, iel, solT4Type);    // global to global mapping between solution node and solution dof
      solT4[i] = (*sol->_Sol[solT4Index])(solT4Dof);      // global extraction and local storage for the solution
      solT4old[i] = (*sol->_SolOld[solT4Index])(solT4Dof);      // global extraction and local storage for the solution
      sysDof[i+nDofsC+nDofsD] = pdeSys->GetSystemDof(solT4Index, solT4PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }
    
    for(unsigned i = 0; i < nDofsT8; i++) {
      unsigned solT8Dof = msh->GetSolutionDof(i, iel, solT8Type);    // global to global mapping between solution node and solution dof
      solT8[i] = (*sol->_Sol[solT8Index])(solT8Dof);      // global extraction and local storage for the solution
      solT8old[i] = (*sol->_SolOld[solT8Index])(solT8Dof);      // global extraction and local storage for the solution
      sysDof[i+nDofsC+nDofsD+nDofsT4] = pdeSys->GetSystemDof(solT8Index, solT8PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }
    
    for(unsigned i = 0; i < nDofsTr; i++) {
      unsigned solTrDof = msh->GetSolutionDof(i, iel, solTrType);    // global to global mapping between solution node and solution dof
      solTr[i] = (*sol->_Sol[solTrIndex])(solTrDof);      // global extraction and local storage for the solution
      solTrold[i] = (*sol->_SolOld[solTrIndex])(solTrDof);      // global extraction and local storage for the solution
      sysDof[i+nDofsC+nDofsD+nDofsT4+nDofsT8] = pdeSys->GetSystemDof(solTrIndex, solTrPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }
    
    for(unsigned i = 0; i < nDofsI2; i++) {
      unsigned solI2Dof = msh->GetSolutionDof(i, iel, solI2Type);    // global to global mapping between solution node and solution dof
      solI2[i] = (*sol->_Sol[solI2Index])(solI2Dof);      // global extraction and local storage for the solution
      solI2old[i] = (*sol->_SolOld[solI2Index])(solI2Dof);      // global extraction and local storage for the solution
      sysDof[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr] = pdeSys->GetSystemDof(solI2Index, solI2PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }
    
    for(unsigned i = 0; i < nDofsI12; i++) {
      unsigned solI12Dof = msh->GetSolutionDof(i, iel, solI12Type);    // global to global mapping between solution node and solution dof
      solI12[i] = (*sol->_Sol[solI12Index])(solI12Dof);      // global extraction and local storage for the solution
      solI12old[i] = (*sol->_SolOld[solI12Index])(solI12Dof);      // global extraction and local storage for the solution
      sysDof[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2] = pdeSys->GetSystemDof(solI12Index, solI12PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }
    
    for(unsigned i = 0; i < nDofsEC; i++) {
      unsigned solECDof = msh->GetSolutionDof(i, iel, solECType);    // global to global mapping between solution node and solution dof
      solEC[i] = (*sol->_Sol[solECIndex])(solECDof);      // global extraction and local storage for the solution
      solECold[i] = (*sol->_SolOld[solECIndex])(solECDof);      // global extraction and local storage for the solution
      sysDof[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2+nDofsI12] = pdeSys->GetSystemDof(solECIndex, solECPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }
    
    for(unsigned i = 0; i < nDofsED; i++) {
      unsigned solEDDof = msh->GetSolutionDof(i, iel, solEDType);    // global to global mapping between solution node and solution dof
      solED[i] = (*sol->_Sol[solEDIndex])(solEDDof);      // global extraction and local storage for the solution
      solEDold[i] = (*sol->_SolOld[solEDIndex])(solEDDof);      // global extraction and local storage for the solution
      sysDof[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2+nDofsI12+nDofsEC] = pdeSys->GetSystemDof(solEDIndex, solEDPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }
    
    for(unsigned i = 0; i < nDofsTalpha; i++) {
      unsigned solTalphaDof = msh->GetSolutionDof(i, iel, solTalphaType);    // global to global mapping between solution node and solution dof
      solTalpha[i] = (*sol->_Sol[solTalphaIndex])(solTalphaDof);      // global extraction and local storage for the solution
      solTalphaold[i] = (*sol->_SolOld[solTalphaIndex])(solTalphaDof);      // global extraction and local storage for the solution
      sysDof[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2+nDofsI12+nDofsEC+nDofsED] = pdeSys->GetSystemDof(solTalphaIndex, solTalphaPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }
    
    for(unsigned i = 0; i < nDofsm1; i++) {
      unsigned solm1Dof = msh->GetSolutionDof(i, iel, solm1Type);    // global to global mapping between solution node and solution dof
      solm1[i] = (*sol->_Sol[solm1Index])(solm1Dof);      // global extraction and local storage for the solution
      solm1old[i] = (*sol->_SolOld[solm1Index])(solm1Dof);      // global extraction and local storage for the solution
      sysDof[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2+nDofsI12+nDofsEC+nDofsED+nDofsTalpha] = pdeSys->GetSystemDof(solm1Index, solm1PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }
    
    for(unsigned i = 0; i < nDofsm2; i++) {
      unsigned solm2Dof = msh->GetSolutionDof(i, iel, solm2Type);    // global to global mapping between solution node and solution dof
      solm2[i] = (*sol->_Sol[solm2Index])(solm2Dof);      // global extraction and local storage for the solution
      solm2old[i] = (*sol->_SolOld[solm2Index])(solm2Dof);      // global extraction and local storage for the solution
      sysDof[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2+nDofsI12+nDofsEC+nDofsED+nDofsTalpha+nDofsm1] = pdeSys->GetSystemDof(solm2Index, solm2PdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }
    
    for(unsigned i = 0; i < nDofsS; i++) {
      unsigned solSDof = msh->GetSolutionDof(i, iel, solSType);    // global to global mapping between solution node and solution dof
      solS[i] = (*sol->_Sol[solSIndex])(solSDof);      // global extraction and local storage for the solution
      solSold[i] = (*sol->_SolOld[solSIndex])(solSDof);      // global extraction and local storage for the solution
      sysDof[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2+nDofsI12+nDofsEC+nDofsED+nDofsTalpha+nDofsm1+nDofsm2] = pdeSys->GetSystemDof(solSIndex, solSPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }
    
    for(unsigned i = 0; i < nDofsu; i++) {
      unsigned soluDof = msh->GetSolutionDof(i, iel, soluType);    // global to global mapping between solution node and solution dof
      solu[i] = (*sol->_Sol[soluIndex])(soluDof);      // global extraction and local storage for the solution
      soluold[i] = (*sol->_SolOld[soluIndex])(soluDof);      // global extraction and local storage for the solution
      sysDof[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2+nDofsI12+nDofsEC+nDofsED+nDofsTalpha+nDofsm1+nDofsm2+nDofsS] = pdeSys->GetSystemDof(soluIndex, soluPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }
    
    for(unsigned i = 0; i < nDofsR; i++) {
      unsigned solRDof = msh->GetSolutionDof(i, iel, solRType);    // global to global mapping between solution node and solution dof
      solR[i] = (*sol->_Sol[solRIndex])(solRDof);      // global extraction and local storage for the solution
      solRold[i] = (*sol->_SolOld[solRIndex])(solRDof);      // global extraction and local storage for the solution
      sysDof[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2+nDofsI12+nDofsEC+nDofsED+nDofsTalpha+nDofsm1+nDofsm2+nDofsS+nDofsu] = pdeSys->GetSystemDof(solRIndex, solRPdeIndex, i, iel);    // global to global mapping between solution node and pdeSys dofs
    }

    // local storage of coordinates
    for(unsigned i = 0; i < nDofsX; i++) {
      unsigned coordXDof  = msh->GetSolutionDof(i, iel, coordXType);    // global to global mapping between coordinates node and coordinate dof
    
//   const double tolerance = 1.e-5;
//   
//   if (face_name == 1) {
//       dirichlet = true;
//         value = 0.;
//   }
//   const double tolerance = 1.e-5;
//   
//   if (face_name == 1) {
//       dirichlet = true;
//         value = 0.;
//   }
      for(unsigned k = 0; k < dim; k++) {
        coordX[k][i] = (*msh->_topology->_Sol[k])(coordXDof);      // global extraction and local storage for the element coordinates
      }
    }

    // start a new recording of all the operations involving adept::adouble variables
    if(assembleMatrix) s.new_recording();

          // *** Face Gauss point loop (boundary Integral) ***
    for (unsigned jface = 0; jface < msh->GetElementFaceNumber (iel); jface++) { //obtain the faces for each elements
      int faceIndex = el->GetBoundaryIndex (iel, jface);
      // look for boundary faces
//printf("%d %d \n", 1000, faceIndex);
      if (faceIndex == 2 || faceIndex == 4) { // faceIndex > 0 indicate that face is exterior. Otherwise, it is interior.
        const unsigned faceGeom = msh->GetElementFaceType (iel, jface);
        unsigned faceDofsT4 = msh->GetElementFaceDofNumber (iel, jface, solT4Type);
        unsigned faceDofsT8 = msh->GetElementFaceDofNumber (iel, jface, solT8Type);
        unsigned faceDofsI12 = msh->GetElementFaceDofNumber (iel, jface, solI12Type);

        vector  < vector  <  double> > faceCoordinatesT4 (dim);   // A matrix holding the face coordinates rowwise.
        vector  < vector  <  double> > faceCoordinatesT8 (dim);   // A matrix holding the face coordinates rowwise.
        vector  < vector  <  double> > faceCoordinatesI12 (dim);   // A matrix holding the face coordinates rowwise.
        
        for (int k = 0; k < dim; k++) {
            faceCoordinatesT4[k].resize (faceDofsT4);
            faceCoordinatesT8[k].resize (faceDofsT8);
            faceCoordinatesI12[k].resize (faceDofsI12);
        }
            
        for (unsigned i = 0; i < faceDofsT4; i++) {
          unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i);   // face-to-element local node mapping.
          for (unsigned k = 0; k < dim; k++) {
            faceCoordinatesT4[k][i] =  coordX[k][inode];
            // We extract the local coordinates on the face from local coordinates on the element.
          }
        }
      
        
        for (unsigned i = 0; i < faceDofsT8; i++) {
          unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i);   // face-to-element local node mapping.
          for (unsigned k = 0; k < dim; k++) {
            faceCoordinatesT8[k][i] =  coordX[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
          }
        }
        
        for (unsigned i = 0; i < faceDofsI12; i++) {
          unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i);   // face-to-element local node mapping.
          for (unsigned k = 0; k < dim; k++) {
            faceCoordinatesI12[k][i] =  coordX[k][inode]; // We extract the local coordinates on the face from local coordinates on the element.
          }
        }
        
        for (unsigned ig = 0; ig  <  msh->_finiteElement[faceGeom][solT4Type]->GetGaussPointNumber(); ig++) {
          // We call the method GetGaussPointNumber from the object finiteElement in the mesh object msh.
          vector < double> normal;
          msh->_finiteElement[faceGeom][solT4Type]->JacobianSur (faceCoordinatesT4, ig, weight, phiT4, phiT4_x, normal);
          msh->_finiteElement[faceGeom][solT8Type]->JacobianSur (faceCoordinatesT8, ig, weight, phiT8, phiT8_x, normal);
          msh->_finiteElement[faceGeom][solI12Type]->JacobianSur (faceCoordinatesI12, ig, weight, phiI12, phiI12_x, normal);

          adept::adouble solT4_gss = 0;
          adept::adouble solT8_gss = 0;
          adept::adouble solI12_gss = 0;
          
          double solT4old_gss = 0;
          double solT8old_gss = 0;
          double solI12old_gss = 0;


          for (unsigned i = 0; i < faceDofsT4; i++) {
            unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i);   // face-to-element local node mapping.
            solT4_gss += phiT4[i] * solT4[inode];
            solT4old_gss += phiT4[i] * solT4old[inode];
          }
          
          for (unsigned i = 0; i < faceDofsT8; i++) {
            unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i);   // face-to-element local node mapping.
            solT8_gss += phiT8[i] * solT8[inode];
            solT8old_gss += phiT8[i] * solT8old[inode];
          }
         
          for (unsigned i = 0; i < faceDofsI12; i++) {
            unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i);   // face-to-element local node mapping.
            solI12_gss += phiI12[i] * solI12[inode];
            solI12old_gss += phiI12[i] * solI12old[inode];
          }
         
          
          for (unsigned i = 0; i < faceDofsT4; i++) {
            unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i);
            aResT4[inode] +=  - alpha * solI12_gss /(solI12_gss + K_I12) * (solT4_gss-T_hat4) * phiT4[i] * weight;    
// printf("%d %d %d %f\n", 1000, faceDofsT4, inode, aResT4[inode]);
          }
          
          for (unsigned i = 0; i < faceDofsT8; i++) {
            unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i);
            aResT8[inode] +=  - alpha * solI12_gss /(solI12_gss + K_I12) * (solT8_gss-T_hat8) * phiT8[i] * weight;    
printf("%d %d %d %e\n", 1000, faceDofsT8, inode, aResT8[inode]);
          }
          
/*         
// aResR[i] += (- (solR_gss - solRold_gss) * phiR[i] / dt - 0.5 * (Temp + TempOld)) * weight;
          
          // *** phi_i loop ***
          double eps = 0.002;
          for (unsigned i = 0; i < faceDofs; i++) {
            unsigned inode = msh->GetLocalFaceVertexIndex (iel, jface, i);
            aRes[inode] +=  phi[i] * eps * (1.0 * solu_gss + 0. * soluOld_gss) * weight;
          }*/
      }
     }
    }
    
    // *** Gauss point loop ***
    for(unsigned ig = 0; ig < msh->_finiteElement[ielGeom][solCType]->GetGaussPointNumber(); ig++) { // guass points for integral
        
      // *** get gauss point weight, test function and test function partial derivatives ***
      msh->_finiteElement[ielGeom][solCType]->Jacobian(coordX, ig, weight, phiC, phiC_x, phiC_xx);
      msh->_finiteElement[ielGeom][solDType]->Jacobian(coordX, ig, weight, phiD, phiD_x, phiD_xx);
      msh->_finiteElement[ielGeom][solT4Type]->Jacobian(coordX, ig, weight, phiT4, phiT4_x, phiT4_xx);
      msh->_finiteElement[ielGeom][solT8Type]->Jacobian(coordX, ig, weight, phiT8, phiT8_x, phiT8_xx);
      msh->_finiteElement[ielGeom][solTrType]->Jacobian(coordX, ig, weight, phiTr, phiTr_x, phiTr_xx);
      msh->_finiteElement[ielGeom][solI2Type]->Jacobian(coordX, ig, weight, phiI2, phiI2_x, phiI2_xx);
      msh->_finiteElement[ielGeom][solI12Type]->Jacobian(coordX, ig, weight, phiI12, phiI12_x, phiI12_xx);
      msh->_finiteElement[ielGeom][solECType]->Jacobian(coordX, ig, weight, phiEC, phiEC_x, phiEC_xx);
      msh->_finiteElement[ielGeom][solEDType]->Jacobian(coordX, ig, weight, phiED, phiED_x, phiED_xx);
      msh->_finiteElement[ielGeom][solTalphaType]->Jacobian(coordX, ig, weight, phiTalpha, phiTalpha_x, phiTalpha_xx);
      msh->_finiteElement[ielGeom][solm1Type]->Jacobian(coordX, ig, weight, phim1, phim1_x, phim1_xx);
      msh->_finiteElement[ielGeom][solm2Type]->Jacobian(coordX, ig, weight, phim2, phim2_x, phim2_xx);
      msh->_finiteElement[ielGeom][solSType]->Jacobian(coordX, ig, weight, phiS, phiS_x, phiS_xx);
      msh->_finiteElement[ielGeom][soluType]->Jacobian(coordX, ig, weight, phiu, phiu_x, phiu_xx);
      msh->_finiteElement[ielGeom][solRType]->Jacobian(coordX, ig, weight, phiR, phiR_x, phiR_xx);
     // msh->_finiteElement[ielGeom][solVType]->Jacobian(coordX, ig, weight, phiV, phiV_x, phiV_xx);

      // evaluate the solution, the solution derivatives and the coordinates in the gauss point
      
      adept::adouble solC_gss = 0.0;
      adept::adouble solD_gss = 0.0;
      adept::adouble solT4_gss = 0.0;
      adept::adouble solT8_gss = 0.0;
      adept::adouble solTr_gss = 0.0;
      adept::adouble solI2_gss = 0.0;
      adept::adouble solI12_gss = 0.0;
      adept::adouble solEC_gss = 0.0;
      adept::adouble solED_gss = 0.0;
      adept::adouble solTalpha_gss = 0.0;
      adept::adouble solm1_gss = 0.0;
      adept::adouble solm2_gss = 0.0;
      adept::adouble solS_gss = 0.0;
      adept::adouble solu_gss = 0.0;
      adept::adouble solR_gss = 0.0;
      
      vector  < double> coordXC_gss;
      vector  < double> coordXD_gss;
      vector  < double> coordXT4_gss;
      vector  < double> coordXT8_gss;
      vector  < double> coordXTr_gss;
      vector  < double> coordXI2_gss;
      vector  < double> coordXI12_gss;
      vector  < double> coordXEC_gss;
      vector  < double> coordXED_gss;
      vector  < double> coordXTalpha_gss;
      vector  < double> coordXm1_gss;
      vector  < double> coordXm2_gss;
      vector  < double> coordXS_gss;
      vector  < double> coordXu_gss;
      vector  < double> coordXR_gss;
      
      coordXC_gss.resize(dim);
      coordXD_gss.resize(dim);
      coordXT4_gss.resize(dim);
      coordXT8_gss.resize(dim);
      coordXTr_gss.resize(dim);
      coordXI2_gss.resize(dim);
      coordXI12_gss.resize(dim);
      coordXEC_gss.resize(dim);
      coordXED_gss.resize(dim);
      coordXTalpha_gss.resize(dim);
      coordXm1_gss.resize(dim);
      coordXm2_gss.resize(dim);
      coordXS_gss.resize(dim);
      coordXu_gss.resize(dim);
      coordXR_gss.resize(dim);
      
      std::fill(coordXC_gss.begin(), coordXC_gss.end(), 0.0);    //set aRes to zero
      std::fill(coordXD_gss.begin(), coordXD_gss.end(), 0.0);    //set aRes to zero
      std::fill(coordXT4_gss.begin(), coordXT4_gss.end(), 0.0);    //set aRes to zero
      std::fill(coordXT8_gss.begin(), coordXT8_gss.end(), 0.0);    //set aRes to zero
      std::fill(coordXTr_gss.begin(), coordXTr_gss.end(), 0.0);    //set aRes to zero
      std::fill(coordXI2_gss.begin(), coordXI2_gss.end(), 0.0);    //set aRes to zero
      std::fill(coordXI12_gss.begin(), coordXI12_gss.end(), 0.0);    //set aRes to zero
      std::fill(coordXEC_gss.begin(), coordXEC_gss.end(), 0.0);    //set aRes to zero
      std::fill(coordXED_gss.begin(), coordXED_gss.end(), 0.0);    //set aRes to zero
      std::fill(coordXTalpha_gss.begin(), coordXTalpha_gss.end(), 0.0);    //set aRes to zero
      std::fill(coordXm1_gss.begin(), coordXm1_gss.end(), 0.0);    //set aRes to zero
      std::fill(coordXm2_gss.begin(), coordXm2_gss.end(), 0.0);    //set aRes to zero
      std::fill(coordXS_gss.begin(),  coordXS_gss.end(), 0.0);    //set aRes to zero
      std::fill(coordXu_gss.begin(),  coordXu_gss.end(), 0.0);    //set aRes to zero
      std::fill(coordXR_gss.begin(),  coordXR_gss.end(), 0.0);    //set aRes to zero
      
      
      double solCold_gss = 0.0;
      double solDold_gss = 0.0;
      double solT4old_gss = 0.0;
      double solT8old_gss = 0.0;
      double solTrold_gss = 0.0;
      double solI2old_gss = 0.0;
      double solI12old_gss = 0.0;
      double solECold_gss = 0.0;
      double solEDold_gss = 0.0;
      double solTalphaold_gss = 0.0;
      double solm1old_gss = 0.0;
      double solm2old_gss = 0.0;
      double solSold_gss = 0.0;
      double soluold_gss = 0.0;
      double solRold_gss = 0.0;
      
      vector < adept::adouble > gradSolC_gss(dim, 0.);
      vector < adept::adouble > gradSolD_gss(dim, 0.);
      vector < adept::adouble > gradSolT4_gss(dim, 0.);
      vector < adept::adouble > gradSolT8_gss(dim, 0.);
      vector < adept::adouble > gradSolTr_gss(dim, 0.);
      vector < adept::adouble > gradSolI2_gss(dim, 0.);
      vector < adept::adouble > gradSolI12_gss(dim, 0.);
      vector < adept::adouble > gradSolEC_gss(dim, 0.);
      vector < adept::adouble > gradSolED_gss(dim, 0.);
      vector < adept::adouble > gradSolTalpha_gss(dim, 0.);
      vector < adept::adouble > gradSolm1_gss(dim, 0.);
      vector < adept::adouble > gradSolm2_gss(dim, 0.);
      vector < adept::adouble > gradSolS_gss(dim, 0.);
      vector < adept::adouble > gradSolu_gss(dim, 0.);
      vector < adept::adouble > gradSolR_gss(dim, 0.);
      
      vector < double > gradSolCold_gss(dim, 0.);
      vector < double > gradSolDold_gss(dim, 0.);
      vector < double > gradSolT4old_gss(dim, 0.);
      vector < double > gradSolT8old_gss(dim, 0.);
      vector < double > gradSolTrold_gss(dim, 0.);
      vector < double > gradSolI2old_gss(dim, 0.);
      vector < double > gradSolI12old_gss(dim, 0.);
      vector < double > gradSolECold_gss(dim, 0.);
      vector < double > gradSolEDold_gss(dim, 0.);
      vector < double > gradSolTalphaold_gss(dim, 0.);
      vector < double > gradSolm1old_gss(dim, 0.);
      vector < double > gradSolm2old_gss(dim, 0.);
      vector < double > gradSolSold_gss(dim, 0.);
      vector < double > gradSoluold_gss(dim, 0.);
      vector < double > gradSolRold_gss(dim, 0.);
      
      
      for(unsigned k = 0; k < dim; k++) {
          for (unsigned i = 0; i<nDofsC;i++){
              coordXC_gss[k] +=coordX[k][i] * phiC[i] ;
          }
          for (unsigned i = 0; i<nDofsD;i++){
              coordXD_gss[k] +=coordX[k][i] * phiD[i] ;
          }
          for (unsigned i = 0; i<nDofsT4;i++){
              coordXT4_gss[k] +=coordX[k][i] * phiT4[i] ;
          }
          for (unsigned i = 0; i<nDofsT8;i++){
              coordXT8_gss[k] +=coordX[k][i] * phiT8[i] ;
          }
          for (unsigned i = 0; i<nDofsTr;i++){
              coordXTr_gss[k] +=coordX[k][i] * phiTr[i]; 
          }
          for (unsigned i = 0; i<nDofsI2;i++){
              coordXI2_gss[k] +=coordX[k][i] * phiI2[i] ;
          }
          for (unsigned i = 0; i<nDofsI12;i++){
              coordXI12_gss[k] +=coordX[k][i] * phiI12[i] ;
          }
          for (unsigned i = 0; i<nDofsEC;i++){
              coordXEC_gss[k] +=coordX[k][i] * phiEC[i] ;
          }
          for (unsigned i = 0; i<nDofsED;i++){
              coordXED_gss[k] +=coordX[k][i] * phiED[i] ;
          }
          for (unsigned i = 0; i<nDofsTalpha;i++){
              coordXTalpha_gss[k] +=coordX[k][i] * phiTalpha[i] ;
          }
          for (unsigned i = 0; i<nDofsm1;i++){
              coordXm1_gss[k] +=coordX[k][i] * phim1[i] ;
          }
          for (unsigned i = 0; i<nDofsm2;i++){
              coordXm2_gss[k] +=coordX[k][i] * phim2[i] ;
          }
          for (unsigned i = 0; i<nDofsS;i++){
              coordXS_gss[k] +=coordX[k][i] * phiS[i] ;
          }
          for (unsigned i = 0; i<nDofsu;i++){
              coordXu_gss[k] +=coordX[k][i] * phiu[i] ;
          }
          for (unsigned i = 0; i<nDofsR;i++){
              coordXR_gss[k] +=coordX[k][i] * phiR[i] ;
          }
      }
      
      unsigned dim1 = dim -1; // for fake dimenion
      
      for(unsigned i = 0; i < nDofsC; i++) {
        solC_gss += phiC[i] * solC[i];
        solCold_gss += phiC[i] * solCold[i];
        for(unsigned j = 0; j < dim; j++) {
          gradSolC_gss[j] += phiC_x[i * dim + j] * solC[i];
          gradSolCold_gss[j] += phiC_x[i * dim + j] * solCold[i];
        }
      }
      
      for(unsigned i = 0; i < nDofsD; i++) {
        solD_gss += phiD[i] * solD[i];
        solDold_gss += phiD[i] * solDold[i];
        for(unsigned j = 0; j < dim; j++) {
          gradSolD_gss[j] += phiD_x[i * dim + j] * solD[i];
          gradSolDold_gss[j] += phiD_x[i * dim + j] * solDold[i];
        }
      }
      
      for(unsigned i = 0; i < nDofsT4; i++) {
        solT4_gss += phiT4[i] * solT4[i];
        solT4old_gss += phiT4[i] * solT4old[i];
        for(unsigned j = 0; j < dim; j++) {
          gradSolT4_gss[j] += phiT4_x[i * dim + j] * solT4[i];
          gradSolT4old_gss[j] += phiT4_x[i * dim + j] * solT4old[i];
        }
      }

     for(unsigned i = 0; i < nDofsT8; i++) {
        solT8_gss += phiT8[i] * solT8[i];
        solT8old_gss += phiT8[i] * solT8old[i];
        for(unsigned j = 0; j < dim; j++) {
          gradSolT8_gss[j] += phiT8_x[i * dim + j] * solT8[i];
          gradSolT8old_gss[j] += phiT8_x[i * dim + j] * solT8old[i];
        }
      }
      
     for(unsigned i = 0; i < nDofsTr; i++) {
        solTr_gss += phiTr[i] * solTr[i];
        solTrold_gss += phiTr[i] * solTrold[i];
        for(unsigned j = 0; j < dim; j++) {
          gradSolTr_gss[j] += phiTr_x[i * dim + j] * solTr[i];
          gradSolTrold_gss[j] += phiTr_x[i * dim + j] * solTrold[i];
        }
      }
      
     for(unsigned i = 0; i < nDofsI2; i++) {
        solI2_gss += phiI2[i] * solI2[i];
        solI2old_gss += phiI2[i] * solI2old[i];
        for(unsigned j = 0; j < dim; j++) {
          gradSolI2_gss[j] += phiI2_x[i * dim + j] * solI2[i];
          gradSolI2old_gss[j] += phiI2_x[i * dim + j] * solI2old[i];
        }
      }
      
     for(unsigned i = 0; i < nDofsI12; i++) {
        solI12_gss += phiI12[i] * solI12[i];
        solI12old_gss += phiI12[i] * solI12old[i];
        for(unsigned j = 0; j < dim; j++) {
          gradSolI12_gss[j] += phiI12_x[i * dim + j] * solI12[i];
          gradSolI12old_gss[j] += phiI12_x[i * dim + j] * solI12old[i];
        }
      }
      
     for(unsigned i = 0; i < nDofsEC; i++) {
        solEC_gss += phiEC[i] * solEC[i];
        solECold_gss += phiEC[i] * solECold[i];
        for(unsigned j = 0; j < dim; j++) {
          gradSolEC_gss[j] += phiEC_x[i * dim + j] * solEC[i];
          gradSolECold_gss[j] += phiEC_x[i * dim + j] * solECold[i];
        }
      }

    for(unsigned i = 0; i < nDofsED; i++) {
        solED_gss += phiED[i] * solED[i];
        solEDold_gss += phiED[i] * solEDold[i];
        for(unsigned j = 0; j < dim; j++) {
          gradSolED_gss[j] += phiED_x[i * dim + j] * solED[i];
          gradSolEDold_gss[j] += phiED_x[i * dim + j] * solEDold[i];
        }
      }
      
    for(unsigned i = 0; i < nDofsTalpha; i++) {
        solTalpha_gss += phiTalpha[i] * solTalpha[i];
        solTalphaold_gss += phiTalpha[i] * solTalphaold[i];
        for(unsigned j = 0; j < dim; j++) {
          gradSolTalpha_gss[j] += phiTalpha_x[i * dim + j] * solTalpha[i];
          gradSolTalphaold_gss[j] += phiTalpha_x[i * dim + j] * solTalphaold[i];
        }
      }
      
    for(unsigned i = 0; i < nDofsm1; i++) {
        solm1_gss += phim1[i] * solm1[i];
        solm1old_gss += phim1[i] * solm1old[i];
        for(unsigned j = 0; j < dim; j++) {
          gradSolm1_gss[j] += phim1_x[i * dim + j] * solm1[i];
          gradSolm1old_gss[j] += phim1_x[i * dim + j] * solm1old[i];
        }
      }
      
    for(unsigned i = 0; i < nDofsm2; i++) {
        solm2_gss += phim2[i] * solm2[i];
        solm2old_gss += phim2[i] * solm2old[i];
        for(unsigned j = 0; j < dim; j++) {
          gradSolm2_gss[j] += phim2_x[i * dim + j] * solm2[i];
          gradSolm2old_gss[j] += phim2_x[i * dim + j] * solm2old[i];
        }
      }
      
    for(unsigned i = 0; i < nDofsS; i++) {
        solS_gss += phiS[i] * solS[i];
        solSold_gss += phiS[i] * solSold[i];
        for(unsigned j = 0; j < dim; j++) {
          gradSolS_gss[j] += phiS_x[i * dim + j] * solS[i];
          gradSolSold_gss[j] += phiS_x[i * dim + j] * solSold[i];
        }
      }
      
    for(unsigned i = 0; i < nDofsu; i++) {
        solu_gss += phiu[i] * solu[i];
        soluold_gss += phiu[i] * soluold[i];
        for(unsigned j = 0; j < dim; j++) {
          gradSolu_gss[j] += phiu_x[i * dim + j] * solu[i];
          gradSoluold_gss[j] += phiu_x[i * dim + j] * soluold[i];
        }
      }
      
    for(unsigned i = 0; i < nDofsR; i++) {
        solR_gss += phiR[i] * solR[i];
        solRold_gss += phiR[i] * solRold[i];
        for(unsigned j = 0; j < dim; j++) {
          gradSolR_gss[j] += phiR_x[i * dim + j] * solR[i];
          gradSolRold_gss[j] += phiR_x[i * dim + j] * solRold[i];
        }
      }

      double dt = mlPdeSys -> GetIntervalTime();
      // *** phiC_i loop ***
      for(unsigned i = 0; i < nDofsC; i++) {
        adept::adouble Temp = 0.;
        adept::adouble TempOld = 0.;
        for(unsigned j = 0; j < dim1; j++) {     
            Temp += 2.0 * solu_gss * solC_gss/coordXC_gss[j] * phiC[i] + solC_gss * gradSolu_gss[j] * phiC[i] + solu_gss * gradSolC_gss[j] *phiC[i];   
            Temp += D_C * gradSolC_gss[j] * phiC_x[i*dim1+j];
            
            TempOld += 2.0 * soluold_gss * solCold_gss/coordXC_gss[j] * phiC[i] + solCold_gss * gradSoluold_gss[j] * phiC[i] 
                    + soluold_gss * gradSolCold_gss[j] *phiC[i];   
            TempOld += D_C * gradSolCold_gss[j] * phiC_x[i*dim1+j];
        }
        Temp += -(lambda_C + lambda_1 * solm1_gss/(solm1_gss + K_m1))*solC_gss*(1.0-solC_gss/C_0)*phiC[i];
        Temp += lambda_T8C * solT8_gss * solC_gss /(1.0 + solEC_gss/K_EC)*phiC[i];
        Temp += (d_C + lambda_S * solS_gss /(solS_gss+K_S))*solC_gss * phiC[i];              
        
        TempOld += -(lambda_C + lambda_1 * solm1_gss/(solm1old_gss + K_m1))*solCold_gss*(1.0-solCold_gss/C_0)*phiC[i];
        TempOld += lambda_T8C * solT8old_gss * solCold_gss /(1.0 + solECold_gss/K_EC)*phiC[i];
        TempOld += (d_C + lambda_S * solSold_gss /(solSold_gss+K_S))*solCold_gss * phiC[i];

        aResC[i] += (- (solC_gss - solCold_gss) * phiC[i] / dt - 0.5 * (Temp + TempOld)) * weight;
      } // end phiC_i loop

      
     for(unsigned i = 0; i < nDofsD; i++) {
        adept::adouble Temp = 0.;
        adept::adouble TempOld = 0.;
        for(unsigned j = 0; j < dim1; j++) {
            Temp    += 2.0 * solu_gss * solD_gss/coordXD_gss[j] * phiD[i]       + solD_gss * gradSolu_gss[j] * phiD[i] 
                    + solu_gss * gradSolD_gss[j] * phiD[i];   
            Temp    += D_D * gradSolD_gss[j] * phiD_x[i*dim1+j];
            
            TempOld += 2.0 * soluold_gss * solDold_gss/coordXD_gss[j] * phiD[i] + solDold_gss * gradSoluold_gss[j] * phiD[i] 
                    + soluold_gss * gradSolDold_gss[j] * phiD[i];   
            TempOld += D_D * gradSolDold_gss[j] * phiD_x[i*dim1+j];
        }
        Temp    += -lambda_D * solC_gss/(solC_gss + K_C) * phiD[i]       + d_D * solD_gss * phiD[i];   
        TempOld += -lambda_D * solCold_gss/(solCold_gss + K_C) * phiD[i] + d_D * solDold_gss * phiD[i];
        
        aResD[i] += (- (solD_gss - solDold_gss) * phiD[i] / dt - 0.5 * (Temp + TempOld)) * weight;
      } // end phiD_i loop
      
      
      for(unsigned i = 0; i < nDofsT4; i++) {
        adept::adouble Temp = 0.;
        adept::adouble TempOld = 0.;
        for(unsigned j = 0; j < dim1; j++) {
            Temp    += 2.0 * solu_gss * solT4_gss/coordXT4_gss[j] * phiT4[i]       + solT4_gss * gradSolu_gss[j] * phiT4[i] 
                    + solu_gss * gradSolT4_gss[j] *phiT4[i];   
            Temp    += D_T4 * gradSolT4_gss[j] * phiT4_x[i*dim1+j];
            
            TempOld += 2.0 * soluold_gss * solT4old_gss/coordXT4_gss[j] * phiT4[i] + solT4old_gss * gradSoluold_gss[j] * phiT4[i] 
                    + soluold_gss * gradSolT4old_gss[j] *phiT4[i];   
            TempOld += D_T4 * gradSolT4old_gss[j] * phiT4_x[i*dim1+j];
        }
        Temp    += - lambda_T4 * T_0 * solI12_gss/(solI12_gss + K_I12)/(1.0 + solTr_gss/K_Tr)/(1.0+solm2_gss/K_m2) * phiT4[i]
                   - lambda_T4I2 * solT4_gss * solI2_gss/(solI2_gss + K_I2) * phiT4[i] 
                   + d_T4 * solT4_gss * phiT4[i];
        TempOld += - lambda_T4 * T_0 * solI12old_gss/(solI12old_gss + K_I12)/(1.0 + solTrold_gss/K_Tr)/(1.0+solm2old_gss/K_m2) * phiT4[i]
                   - lambda_T4I2 * solT4old_gss * solI2old_gss/(solI2old_gss + K_I2) * phiT4[i]
                   + d_T4 * solT4old_gss * phiT4[i];
        
        aResT4[i] += (- (solT4_gss - solT4old_gss) * phiT4[i] / dt - 0.5 * (Temp + TempOld)) * weight;
      } // end phiT4_i loop
      
      for(unsigned i = 0; i < nDofsT8; i++) {
        adept::adouble Temp = 0.;
        adept::adouble TempOld = 0.;
        for(unsigned j = 0; j < dim1; j++) {
            Temp    += 2.0 * solu_gss * solT8_gss/coordXT8_gss[j] * phiT8[i] + solT8_gss * gradSolu_gss[j] * phiT8[i] 
                    + solu_gss * gradSolT8_gss[j] * phiT8[i];   
            Temp    += D_T8 * gradSolT8_gss[j] * phiT8_x[i*dim1+j];
            
            TempOld += 2.0 * soluold_gss * solT8old_gss/coordXT8_gss[j] * phiT8[i] + solT8old_gss * gradSoluold_gss[j] * phiT8[i] 
                    + soluold_gss * gradSolT8old_gss[j] *phiT8[i];   
            TempOld += D_T8 * gradSolT8old_gss[j] * phiT8_x[i*dim1+j];
        }
        Temp    += - lambda_T8 * T_80 * solI12_gss/(solI12_gss + K_I12)/(1.0 + solTr_gss/K_Tr)/(1.0+solm2_gss/K_m2) * phiT8[i]
                   - lambda_T8I2 * solT8_gss * solI2_gss/(solI2_gss + K_I2) * phiT8[i]
                   + d_T8 * solT8_gss * phiT8[i];
        TempOld += - lambda_T8 * T_80 * solI12old_gss/(solI12old_gss + K_I12)/(1.0 + solTrold_gss/K_Tr)/(1.0+solm2old_gss/K_m2) * phiT8[i]
                   - lambda_T8I2 * solT8old_gss * solI2old_gss/(solI2old_gss + K_I2) * phiT8[i]
                   + d_T8 * solT8old_gss * phiT8[i] ;
        
        aResT8[i] += (- (solT8_gss - solT8old_gss) * phiT8[i] / dt - 0.5 * (Temp + TempOld)) * weight;
      } // end phiT8_i loop
      
      
      for(unsigned i = 0; i < nDofsTr; i++) {
        adept::adouble Temp = 0.;
        adept::adouble TempOld = 0.;
        for(unsigned j = 0; j < dim1; j++) {
            Temp    += 2.0 * solu_gss * solTr_gss/coordXTr_gss[j] * phiTr[i] + solTr_gss * gradSolu_gss[j] * phiTr[i] 
                    +  solu_gss * gradSolTr_gss[j] *phiTr[i];   
            Temp    += D_Tr * gradSolTr_gss[j] * phiTr_x[i*dim1+j];
            
            TempOld += 2.0 * soluold_gss * solTrold_gss/coordXTr_gss[j] * phiTr[i] + solTrold_gss * gradSoluold_gss[j] * phiTr[i] 
                    +  soluold_gss * gradSolTrold_gss[j] *phiTr[i];   
            TempOld += D_Tr * gradSolTrold_gss[j] * phiTr_x[i*dim1+j];
        }
        Temp    += - lambda_Tr * T_0 * solEC_gss/(solEC_gss + K_EC)/(1.0 + solTalpha_gss/K_Talpha) * phiTr[i]          + d_Tr * solTr_gss * phiTr[i]; 
        TempOld += - lambda_Tr * T_0 * solECold_gss/(solECold_gss + K_EC)/(1.0 + solTalphaold_gss/K_Talpha) * phiTr[i] + d_Tr * solTrold_gss * phiTr[i];
       
        aResTr[i] += (- (solTr_gss - solTrold_gss) * phiTr[i] / dt - 0.5 * (Temp + TempOld)) * weight;
      } // end phiTr_i loop
      
      for(unsigned i = 0; i < nDofsI2; i++) {
        adept::adouble Temp = 0.;
        adept::adouble TempOld = 0.;
        for(unsigned j = 0; j < dim1; j++) {
            Temp    += D_I2 * gradSolI2_gss[j] * phiI2_x[i*dim1+j];

            TempOld += D_I2 * gradSolI2old_gss[j] * phiI2_x[i*dim1+j];
        }
        Temp    += -lambda_I2T4 * solT4_gss * phiI2[i]    + d_I2 * solI2_gss * phiI2[i];
        TempOld += -lambda_I2T4 * solT4old_gss * phiI2[i] + d_I2 * solI2old_gss * phiI2[i];
        
        aResI2[i] += (- (solI2_gss - solI2old_gss) * phiI2[i] / dt - 0.5 * (Temp + TempOld)) * weight;
      } // end phiI2_i loop
      
      for(unsigned i = 0; i < nDofsI12; i++) {
        adept::adouble Temp = 0.;
        adept::adouble TempOld = 0.;
        for(unsigned j = 0; j < dim1; j++) {
            Temp    += D_I12 * gradSolI12_gss[j]    * phiI12_x[i*dim1+j];

            TempOld += D_I12 * gradSolI12old_gss[j] * phiI12_x[i*dim1+j];
        }
        Temp    += -lambda_I12D * solD_gss    * phiI12[i] + d_I12 * solI12_gss    * phiI12[i];
        TempOld += -lambda_I12D * solDold_gss * phiI12[i] + d_I12 * solI12old_gss * phiI12[i];
        
        aResI12[i] += (- (solI12_gss - solI12old_gss) * phiI12[i] / dt - 0.5 * (Temp + TempOld)) * weight;
      } // end phiI12_i loop
      
     for(unsigned i = 0; i < nDofsEC; i++) {
        adept::adouble Temp = 0.;
        adept::adouble TempOld = 0.;
        for(unsigned j = 0; j < dim1; j++) {
            Temp    += D_EC * gradSolEC_gss[j]    * phiEC_x[i*dim1+j];

            TempOld += D_EC * gradSolECold_gss[j] * phiEC_x[i*dim1+j];
        }
        Temp    += -lambda_EC * solC_gss * phiEC[i];
        Temp    += d_EC * solEC_gss * (solT4_gss/(solT4_gss + K_T4)         + solT8_gss/(solT8_gss + K_T8)        + solTr_gss/(solTr_gss + K_Tr)
                + solD_gss/(solD_gss + K_D)       + solC_gss/(solC_gss + K_C))      * phiEC[i];
        
        TempOld += -lambda_EC * solCold_gss * phiEC[i];
        TempOld += d_EC * solECold_gss * (solT4old_gss/(solT4old_gss + K_T4) + solT8old_gss/(solT8old_gss + K_T8) + solTrold_gss/(solTrold_gss + K_Tr)
                + solDold_gss/(solDold_gss + K_D) + solCold_gss/(solCold_gss + K_C)) * phiEC[i];
        
        aResEC[i] += (- (solEC_gss - solECold_gss) * phiEC[i] / dt - 0.5 * (Temp + TempOld)) * weight;
      } // end phiEC_i loop
      
    for(unsigned i = 0; i < nDofsED; i++) {
        adept::adouble Temp = 0.;
        adept::adouble TempOld = 0.;
        for(unsigned j = 0; j < dim1; j++) {
            Temp    += D_ED * gradSolED_gss[j] * phiED_x[i*dim1+j];

            TempOld += D_ED * gradSolECold_gss[j] * phiED_x[i*dim1+j];
        }
        Temp    += -lambda_ED * solD_gss * phiED[i];
        Temp    += d_ED * solED_gss * (solTr_gss/(solTr_gss + K_Tr)          + solC_gss/(solC_gss + K_C)) * phiED[i];
        
        TempOld += -lambda_ED * solDold_gss * phiED[i];
        TempOld += d_ED * solEDold_gss * (solTrold_gss/(solTrold_gss + K_Tr) + solCold_gss/(solCold_gss + K_C)) * phiED[i];
        
        aResED[i] += (- (solED_gss - solEDold_gss) * phiED[i] / dt - 0.5 * (Temp + TempOld)) * weight;
      } // end phiED_i loop
      
      for(unsigned i = 0; i < nDofsTalpha; i++) {
        adept::adouble Temp = 0.;
        adept::adouble TempOld = 0.;
        for(unsigned j = 0; j < dim1; j++) {
            Temp    += D_Talpha * gradSolTalpha_gss[j] * phiTalpha_x[i*dim1+j];

            TempOld += D_Talpha * gradSolTalphaold_gss[j] * phiTalpha_x[i*dim1+j];
        }
        Temp    += -lambda_TalphaD  * solD_gss  * phiTalpha[i];
        Temp    += -lambda_TalphaED * solED_gss * solTr_gss/(solTr_gss + K_Tr) * phiTalpha[i]; 
        Temp    += d_Talpha * solTalpha_gss * phiTalpha[i];
        
        TempOld += -lambda_TalphaD * solDold_gss * phiTalpha[i];
        TempOld += -lambda_TalphaED * solEDold_gss * solTrold_gss/(solTrold_gss + K_Tr) * phiTalpha[i]; 
        TempOld += d_Talpha * solTalphaold_gss * phiTalpha[i]  ;

        aResTalpha[i] += (- (solTalpha_gss - solTalphaold_gss) * phiTalpha[i] / dt - 0.5 * (Temp + TempOld)) * weight;
      } // end phiTalpha_i loop
      
      for(unsigned i = 0; i < nDofsm1; i++) {
        adept::adouble Temp = 0.;
        adept::adouble TempOld = 0.;
        for(unsigned j = 0; j < dim1; j++) {
            Temp    += D_m1 * gradSolm1_gss[j] * phim1_x[i*dim1+j];

            TempOld += D_m1 * gradSolm1old_gss[j] * phim1_x[i*dim1+j];
        }
        Temp    += -lambda_m1EC * solEC_gss    * (solC_gss/(solC_gss + K_C)) * phim1[i]       + d_m1 * solm1_gss * phim1[i];
        
        TempOld += -lambda_m1EC * solECold_gss * (solCold_gss/(solCold_gss + K_C)) * phim1[i] + d_m1 * solm1old_gss * phim1[i];

        aResm1[i] += (- (solm1_gss - solm1old_gss) * phim1[i] / dt - 0.5 * (Temp + TempOld)) * weight;
      } // end phim1 loop
      
      for(unsigned i = 0; i < nDofsm2; i++) {
        adept::adouble Temp = 0.;
        adept::adouble TempOld = 0.;
        for(unsigned j = 0; j < dim1; j++) {
            Temp    += D_m2 * gradSolm2_gss[j] * phim2_x[i*dim1+j];

            TempOld += D_m2 * gradSolm2old_gss[j] * phim2_x[i*dim1+j];
        }
        Temp      += -lambda_m2EC * solEC_gss * (solD_gss/(solD_gss + K_D)) * phim2[i]          + d_m2 * solm2_gss * phim2[i];
        
        TempOld   += -lambda_m2EC * solECold_gss * (solDold_gss/(solDold_gss + K_D)) * phim2[i] + d_m2 * solm2old_gss * phim2[i];

        aResm2[i] += (- (solm2_gss - solm2old_gss) * phim2[i] / dt - 0.5 * (Temp + TempOld)) * weight;
      } // end phim2 loop
      
      for(unsigned i = 0; i < nDofsS; i++) {
        adept::adouble Temp = 0.;
        adept::adouble TempOld = 0.;
        for(unsigned j = 0; j < dim1; j++) {
            Temp    += D_S * gradSolS_gss[j] * phiS_x[i*dim1+j];

            TempOld += D_S * gradSolSold_gss[j] * phiS_x[i*dim1+j];
        }
        Temp    += - A_S * phiS[i] -lambda_SED * solED_gss * solC_gss/(solC_gss + K_C) * phiS[i]          + d_S * solS_gss * phiS[i]; 
        
        TempOld += - A_S * phiS[i] -lambda_SED * solEDold_gss * solCold_gss/(solCold_gss + K_C) * phiS[i] + d_S * solSold_gss * phiS[i];

        aResS[i] += (- (solS_gss - solSold_gss) * phiS[i] / dt - 0.5 * (Temp + TempOld)) * weight;
      } // end phiS_i loop
      
      for(unsigned i = 0; i < nDofsu; i++) {
        adept::adouble Temp = 0.;
        for(unsigned j = 0; j < dim1; j++) {
            Temp += 0.403 * (2.0 * solu_gss/coordXu_gss[j] + gradSolu_gss[j]) * phiu[i]; 
        }
        Temp += -(lambda_C + lambda_1 * solm1_gss/(solm1_gss + K_m1))*solC_gss*(1.0-solC_gss/C_0)*phiu[i];
        Temp +=  lambda_T8C * solT8_gss * solC_gss /(1.0 + solEC_gss/K_EC)*phiu[i];
        Temp +=  (d_C + lambda_S * solS_gss /(solS_gss+K_S))*solC_gss * phiu[i];
        Temp += - lambda_D * solC_gss/(solC_gss + K_C) * phiu[i] + d_D * solD_gss * phiu[i];   
        Temp += - lambda_T4 * T_0 * solI12_gss/(solI12_gss + K_I12)/(1.0 + solTr_gss/K_Tr)/(1.0+solm2_gss/K_m2) * phiu[i]
                - lambda_T4I2 * solT4_gss * solI2_gss/(solI2_gss + K_I2) * phiu[i] 
                + d_T4 * solT4_gss * phiu[i];
        Temp += - lambda_T8 * T_80 * solI12_gss/(solI12_gss + K_I12)/(1.0 + solTr_gss/K_Tr)/(1.0+solm2_gss/K_m2) * phiu[i]
                - lambda_T8I2 * solT8_gss * solI2_gss/(solI2_gss + K_I2) * phiu[i]
                + d_T8 * solT8_gss * phiu[i];
        Temp += - lambda_Tr * T_0 * solEC_gss/(solEC_gss + K_EC)/(1.0 + solTalpha_gss/K_Talpha) * phiu[i] + d_Tr * solTr_gss * phiu[i]; 
        
        aResu[i] += -Temp * weight;
      } // end phiu_i loop
          
    for(unsigned i = 0; i < nDofsR; i++) {
        adept::adouble Temp = 0.;
        adept::adouble TempOld = 0.;
        Temp     += -solu_gss * phiR[i]; 
        TempOld  += -soluold_gss * phiR[i]; 
        aResR[i] += (- (solR_gss - solRold_gss) * phiR[i] / dt - 0.5 * (Temp + TempOld)) * weight;
      } // end phiS_i loo
    } // end gauss point loop

    //--------------------------------------------------------------------------------------------------------
    // Add the local Matrix/Vector into the global Matrix/Vector

    //copy the value of the adept::adoube aRes in double Res and store them in RES
    Res.resize(nDofsTot);    //resize

    for(int i = 0; i < nDofsC; i++) Res[i] = -aResC[i].value();
    
    for(int i = 0; i < nDofsD; i++) Res[i+nDofsC] = -aResD[i].value();
    
    for(int i = 0; i < nDofsT4; i++) Res[i+nDofsC+nDofsD] = -aResT4[i].value();
    
    for(int i = 0; i < nDofsT8; i++) Res[i+nDofsC+nDofsD+nDofsT4] = -aResT8[i].value();

    for(int i = 0; i < nDofsTr; i++) Res[i+nDofsC+nDofsD+nDofsT4+nDofsT8] = -aResTr[i].value();
    
    for(int i = 0; i < nDofsI2; i++) Res[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr] = -aResI2[i].value();
                                                      
    for(int i = 0; i < nDofsI12; i++) Res[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2] = -aResI12[i].value();
    
    for(int i = 0; i < nDofsEC; i++) Res[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2+nDofsI12] = -aResEC[i].value();
    
    for(int i = 0; i < nDofsED; i++) Res[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2+nDofsI12+nDofsEC] = -aResED[i].value();
    
    for(int i = 0; i < nDofsTalpha; i++) Res[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2+nDofsI12+nDofsEC+nDofsED] = -aResTalpha[i].value();
                                                      
    for(int i = 0; i < nDofsm1; i++) Res[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2+nDofsI12+nDofsEC+nDofsED+nDofsTalpha] = -aResm1[i].value();
    
    for(int i = 0; i < nDofsm2; i++) Res[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2+nDofsI12+nDofsEC+nDofsED+nDofsTalpha+nDofsm1] = -aResm2[i].value();
    
    for(int i = 0; i < nDofsS; i++) Res[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2+nDofsI12+nDofsEC+nDofsED+nDofsTalpha+nDofsm1+nDofsm2] = -aResS[i].value();
    
    for(int i = 0; i < nDofsu; i++) Res[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2+nDofsI12+nDofsEC+nDofsED+nDofsTalpha+nDofsm1+nDofsm2+nDofsS] = -aResu[i].value();
    
    for(int i = 0; i < nDofsR; i++) Res[i+nDofsC+nDofsD+nDofsT4+nDofsT8+nDofsTr+nDofsI2+nDofsI12+nDofsEC+nDofsED+nDofsTalpha+nDofsm1+nDofsm2+nDofsS+nDofsu] = -aResR[i].value();
                                                      
    RES->add_vector_blocked(Res, sysDof);

    //Extarct and store the Jacobian
    if(assembleMatrix) {
      Jac.resize(nDofsTot * nDofsTot);  
      // define the dependent variables
      s.dependent(&aResC[0], nDofsC);
      s.dependent(&aResD[0], nDofsD);
      s.dependent(&aResT4[0], nDofsT4);
      s.dependent(&aResT8[0], nDofsT8);
      s.dependent(&aResTr[0], nDofsTr);
      s.dependent(&aResI2[0], nDofsI2);
      s.dependent(&aResI12[0], nDofsI12);
      s.dependent(&aResEC[0], nDofsEC);
      s.dependent(&aResED[0], nDofsED);
      s.dependent(&aResTalpha[0], nDofsTalpha);
      s.dependent(&aResm1[0], nDofsm1);
      s.dependent(&aResm2[0], nDofsm2);
      s.dependent(&aResS[0], nDofsS);
      s.dependent(&aResu[0], nDofsu);
      s.dependent(&aResR[0], nDofsR);
      
      s.independent(&solC[0], nDofsC);
      s.independent(&solD[0], nDofsD);
      s.independent(&solT4[0], nDofsT4);
      s.independent(&solT8[0], nDofsT8);
      s.independent(&solTr[0], nDofsTr);
      s.independent(&solI2[0], nDofsI2);
      s.independent(&solI12[0], nDofsI12);
      s.independent(&solEC[0], nDofsEC);
      s.independent(&solED[0], nDofsED);
      s.independent(&solTalpha[0], nDofsTalpha);
      s.independent(&solm1[0], nDofsm1);
      s.independent(&solm2[0], nDofsm2);
      s.independent(&solS[0], nDofsS);
      s.independent(&solu[0], nDofsu);
      s.independent(&solR[0], nDofsR);

      // get the and store jacobian matrix (row-major)
      s.jacobian(&Jac[0] , true);
      KK->add_matrix_blocked(Jac, sysDof, sysDof);

      s.clear_independents();
      s.clear_dependents();
    }
  } //end element loop for each process

  RES->close();

  if(assembleMatrix) {
    KK->close();
  }

  // ***************** END ASSEMBLY *******************
}

