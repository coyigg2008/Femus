#include "MultiLevelProblem.hpp"
#include "MultiLevelMesh.hpp"
#include "TransientSystem.hpp"
#include "NumericVector.hpp"
#include "Fluid.hpp"
#include "Parameter.hpp"
#include "FemTTUInit.hpp"
#include "SparseMatrix.hpp"
#include "VTKWriter.hpp"
#include "GMVWriter.hpp"
#include "NonLinearImplicitSystem.hpp"
#include "SolvertypeEnum.hpp"
#include "FElemTypeEnum.hpp"

using std::cout;
using std::endl;

using namespace femus;

void AssembleMatrixResNS(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix);
void AssembleMatrixResNS_old(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix);
void AssembleMatrixResT(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix);


double InitVariableU(const double &x, const double &y, const double &z);


bool SetBoundaryConditionTurek(const double &x, const double &y, const double &z,const char name[], 
			       double &value, const int FaceName, const double time);

bool SetBoundaryConditionCavityFlow(const double &x, const double &y, const double &z,const char name[], 
				    double &value, const int FaceName, const double time);

bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber,const int &level);

int main(int argc,char **args) {

  bool Vanka=0, Gmres=0, Asm=0;
  if(argc >= 2) {
    if( !strcmp("vanka",args[1])) 	Vanka=1;
    else if( !strcmp("gmres",args[1])) 	Gmres=1;
    else if( !strcmp("asm",args[1])) 	Asm=1;
    
    if(Vanka+Gmres+Asm==0) {
      cout << "wrong input arguments!" << endl;
      exit(0);
    }
  }
  else {
    cout << "No input argument set default smoother = Gmres" << endl;
    Gmres=1;
  }
  
  /// Init Petsc-MPI communicator
  FemTTUInit mpinit(argc,args,MPI_COMM_WORLD);
  
  /// INIT MESH =================================  
  
  unsigned short nm,nr;
  nm=4;
  std::cout<<"MULTIGRID levels: "<< nm << endl;

  nr=0;
  std::cout<<"MAX_REFINEMENT levels: " << nr << endl<< endl;
  
  int tmp=nm;  nm+=nr;  nr=tmp;
  
  char *infile = new char [50];
 
  sprintf(infile,"./input/box10x10.neu");
  
  //Adimensional quantity (Lref,Uref)
  double Lref = 1.;
  double Uref = 1.;
  
  //Steadystate NonLinearMultiLevelProblem  
  //MultiLevelMesh ml_msh(nm,nr,infile,"seventh",Lref,SetRefinementFlag); 
  MultiLevelMesh ml_msh;
  ml_msh.ReadCoarseMesh(infile,"seventh",Lref);
  ml_msh.RefineMesh(nm,nr,SetRefinementFlag);
  
  ml_msh.EraseCoarseLevels(nm-1);
  
  MultiLevelSolution ml_sol(&ml_msh);
  
  // generate solution vector
//   ml_sol.AddSolution("T",LAGRANGE,SECOND);
   ml_sol.AddSolution("U",LAGRANGE,SECOND);
   ml_sol.AddSolution("V",LAGRANGE,SECOND);
//   // the pressure variable should be the last for the Schur decomposition
//   ml_sol.AddSolution("P",DISCONTINOUS_POLYNOMIAL,FIRST);
  
  ml_sol.AddSolution("T",LAGRANGE,FIRST);
//  ml_sol.AddSolution("U",LAGRANGE,FIRST);
//  ml_sol.AddSolution("V",LAGRANGE,FIRST);
  // the pressure variable should be the last for the Schur decomposition
  ml_sol.AddSolution("P",LAGRANGE,FIRST);
  ml_sol.AssociatePropertyToSolution("P","Pressure");
 
  //Initialize (update Init(...) function)
  ml_sol.Initialize("U",InitVariableU);
  ml_sol.Initialize("V");
  ml_sol.Initialize("P");
  ml_sol.Initialize("T");
  
  //Set Boundary (update Dirichlet(...) function)
  ml_sol.AttachSetBoundaryConditionFunction(SetBoundaryConditionCavityFlow);
  ml_sol.GenerateBdc("U");
  ml_sol.GenerateBdc("V");
  ml_sol.GenerateBdc("P");
  ml_sol.GenerateBdc("T");
  
  MultiLevelProblem ml_prob(&ml_msh,&ml_sol);
  
  // add fluid material
  Parameter parameter(Lref,Uref);
  
  // Generate fluid Object (Adimensional quantities,viscosity,density,fluid-model)
  Fluid fluid(parameter,0.001,1,"Newtonian",0.001,1.);
  cout << "Fluid properties: " << endl;
  cout << fluid << endl;
  
  ml_prob.parameters.set<Fluid>("Fluid") = fluid;
   
  
  //BEGIN Navier-Stokes Multilevel Problem
  std::cout << std::endl;
  std::cout << " *********** Navier-Stokes ************  " << std::endl;
    
  NonLinearImplicitSystem & system1 = ml_prob.add_system<NonLinearImplicitSystem> ("Navier-Stokes");
  system1.AddSolutionToSytemPDE("U");
  system1.AddSolutionToSytemPDE("V");
  system1.AddSolutionToSytemPDE("P");
  
  // Set MG Options
  system1.AttachAssembleFunction(AssembleMatrixResNS);  
  system1.SetMaxNumberOfNonLinearIterations(20);
  system1.SetMaxNumberOfLinearIterations(2);
  system1.SetAbsoluteConvergenceTolerance(1.e-10);
  system1.SetNonLinearConvergenceTolerance(1.e-8);
  system1.SetMgType(F_CYCLE);
  system1.SetNumberPreSmoothingStep(1);
  system1.SetNumberPostSmoothingStep(1);
      
  //Set Smoother Options
  if(Gmres) 		system1.SetMgSmoother(GMRES_SMOOTHER);
  else if(Asm) 		system1.SetMgSmoother(ASM_SMOOTHER);
  else if(Vanka)	system1.SetMgSmoother(VANKA_SMOOTHER);
  
  system1.init();
  //common smoother options
//   system1.AddStabilization(true);
  system1.SetSolverFineGrids(GMRES);
  system1.SetPreconditionerFineGrids(ILU_PRECOND); 
  system1.SetTolerances(1.e-12,1.e-20,1.e+50,4);
 
  system1.ClearVariablesToBeSolved();
  //system1.AddVariableToBeSolved("All");
  system1.AddVariableToBeSolved("U");
  system1.AddVariableToBeSolved("V");
  system1.AddVariableToBeSolved("P");
  //for Vanka and ASM smoothers
  system1.SetNumberOfSchurVariables(0);
  system1.SetElementBlockNumber(4);   
  //system1.SetElementBlockNumber("All",1);     
  //for Gmres smoother
  system1.SetDirichletBCsHandling(PENALTY); 
  //system1.SetDirichletBCsHandling(ELIMINATION); 
   
  // Solve Navier-Stokes system
  ml_prob.get_system("Navier-Stokes").solve();
  //END Navier-Stokes Multilevel Problem
  
  
  //BEGIN Temperature MultiLevel Problem
  std::cout << std::endl;
  std::cout << " *********** Temperature ************* " << std::endl;
    
  LinearImplicitSystem & system2 = ml_prob.add_system<LinearImplicitSystem> ("Temperature");
  system2.AddSolutionToSytemPDE("T");
  
  
  // Set MG Options
  system2.AttachAssembleFunction(AssembleMatrixResT);
  system2.SetMaxNumberOfLinearIterations(6);
  system2.SetAbsoluteConvergenceTolerance(1.e-9);  
  system2.SetMgType(V_CYCLE);
  system2.SetNumberPreSmoothingStep(1);
  system2.SetNumberPostSmoothingStep(1);
   
  //Set Smoother Options
  if(Gmres) 		system2.SetMgSmoother(GMRES_SMOOTHER);
  else if(Asm) 		system2.SetMgSmoother(ASM_SMOOTHER);
  else if(Vanka)	system2.SetMgSmoother(VANKA_SMOOTHER);
  
  system2.init(); 
  //common smoother option
  system2.SetSolverFineGrids(GMRES); 
  system2.SetTolerances(1.e-12,1.e-20,1.e+50,4);
  system2.SetPreconditionerFineGrids(ILU_PRECOND);
  //for Vanka and ASM smoothers
  system2.ClearVariablesToBeSolved();
  system2.AddVariableToBeSolved("All");
  system2.SetNumberOfSchurVariables(0);
  system2.SetElementBlockNumber(4);                
  //for Gmres smoother
  system2.SetDirichletBCsHandling(PENALTY); 
  //system2.SetDirichletBCsHandling(ELIMINATION); 
  
  
  // Solve Temperature system
  ml_prob.get_system("Temperature").solve();
  //END Temperature Multilevel Problem
    
  /// Print all solutions
  std::vector<std::string> print_vars;
  print_vars.push_back("U");
  print_vars.push_back("V");
  print_vars.push_back("P");
  print_vars.push_back("T");
  
     
  VTKWriter vtkio(ml_sol);
  vtkio.write_system_solutions("biquadratic",print_vars);
  
//   XDMFWriter xdmfio(ml_prob);
//   xdmfio.write_system_solutions("biquadratic",print_vars);
  
  GMVWriter gmvio(ml_sol);
  gmvio.write_system_solutions("linear",print_vars);
  
  //Destroy all the new systems
  ml_prob.clear();
  
  delete [] infile;
  return 0;
}

//-----------------------------------------------------------------------------------------------------------------

bool SetRefinementFlag(const double &x, const double &y, const double &z, const int &ElemGroupNumber, const int &level) {
  bool refine=0;
  // refinemenet based on Elemen Group Number
  if(ElemGroupNumber==5 ) {
    refine=1;
  }
  if(ElemGroupNumber==6 && level<2) {
    refine=1;
  }
  if(ElemGroupNumber==7 ) {
    refine=0;
  }

  return refine;
}

//--------------------------------------------------------------------------------------------------------------

double InitVariableU(const double &x, const double &y, const double &z) { 
   double um = 0.2;
   double  value=1.5*um*(4.0/(0.1681))*y*(0.41-y); 
   return value;
}

//-------------------------------------------------------------------------------------------------------------------

bool SetBoundaryConditionTurek(const double &x, const double &y, const double &z,const char name[], 
			       double &value, const int FaceName, const double time){
  bool test=1; //Dirichlet
  value=0.;
  //   cout << "Time bdc : " <<  time << endl;
  if(!strcmp(name,"U")) {
    if(1==FaceName){   //inflow
      test=1;
      double um = 0.2; // U/Uref
      value=1.5*0.2*(4.0/(0.1681))*y*(0.41-y);
    }  
    else if(2==FaceName ){  //outflow
      test=0;
      //    test=1;
      value=0.;
    }
    else if(3==FaceName ){  // no-slip fluid wall
      test=1;
      value=0.;	
    }
    else if(4==FaceName ){  // no-slip solid wall
      test=1;
      value=0.;
    }
  }  
  else if(!strcmp(name,"V")){
    if(1==FaceName){            //inflow
      test=1;
      value=0.;
    }  
    else if(2==FaceName ){      //outflow
      test=0;
      //    test=1;
      value=0.;
    }
    else if(3==FaceName ){      // no-slip fluid wall
      test=1;
      value=0;
    }
    else if(4==FaceName ){      // no-slip solid wall
      test=1;
      value=0.;
    }
  }
  else if(!strcmp(name,"W")){
    if(1==FaceName){
      test=1;
      value=0.;
    }  
    else if(2==FaceName ){  
      test=1;
      value=0.;
    }
    else if(3==FaceName ){  
      test=1;
      value=0.;
    }
    else if(4==FaceName ){  
      test=1;
      value=0.;
    }
  }
  else if(!strcmp(name,"P")){
    if(1==FaceName){
      test=0;
      value=0.;
    }  
    else if(2==FaceName ){  
      test=0;
      value=0.;
    }
    else if(3==FaceName ){  
      test=0;
      value=0.;
    }
    else if(4==FaceName ){  
      test=0;
      value=0.;
    }
  }
  else if(!strcmp(name,"T")) {
    if(1==FaceName){   //inflow
      test=1;
      value=1;
    }  
    else if(2==FaceName ){  //outflow
      test=0;
      value=0.;
    }
    else if(3==FaceName ){  // no-slip fluid wall
      test=0;
      value=0.;	
    }
    else if(4==FaceName ){  // no-slip solid wall
      test=1;
      value=5.;
    }
  }  
  
  return test;
}


bool SetBoundaryConditionCavityFlow(const double& x, const double& y, const double& z, const char name[], double& value, const int FaceName, const double time){
  bool test=1; //Dirichlet
  value=0.;
  if(!strcmp(name,"V")){
    if(1==FaceName){            //inflow
      test=1;
      if(y<0.5 && y>-0.5) value=1.;//4*(0.5-y)*(y+0.5);
    }
  }  
  return test;
}













// //------------------------------------------------------------------------------------------------------------

#include "adept.h"

static unsigned counter=0;

void AssembleMatrixResNS(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix){
 
       
    clock_t AssemblyTime=0;
    clock_t start_time, end_time;
  
    static adept::Stack adeptStack; 
    
    //pointers and references
    MultiLevelSolution*	 ml_sol	                      = ml_prob._ml_sol;
    Solution*	 mysolution  	                      = ml_sol->GetSolutionLevel(level);
    NonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<NonLinearImplicitSystem>("Navier-Stokes");
    LinearEquationSolver*  myLinEqSolver	     = my_nnlin_impl_sys._LinSolver[level];   
    
    mesh		*mymsh		=  ml_prob._ml_msh->GetLevel(level);
    elem		*myel		=  mymsh->el;
    SparseMatrix	*myKK		=  myLinEqSolver->_KK;
    NumericVector 	*myRES		=  myLinEqSolver->_RES;
        
    const unsigned dim = mymsh->GetDimension();
    const unsigned nabla_dim = 3*(dim-1);
    const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
  
    // local objects
    vector<adept::adouble> SolVAR(dim+1);
    vector<vector<adept::adouble> > GradSolVAR(dim+1);
    vector<vector<adept::adouble> > NablaSolVAR(dim+1);
    
    for(int i=0;i<dim+1;i++){
      GradSolVAR[i].resize(dim);
      NablaSolVAR[i].resize(nabla_dim);
    }
        
    vector <double > phi;
    vector <adept::adouble> gradphi;
    vector <adept::adouble> nablaphi;
    adept::adouble Weight;
    
    phi.reserve(max_size);
    gradphi.reserve(max_size*dim);
    nablaphi.reserve(max_size*nabla_dim);
        
    vector <double > phi1;
    vector <adept::adouble> gradphi1;
    vector <adept::adouble> nablaphi1;
    adept::adouble Weight1;
    
    phi1.reserve(max_size);
    gradphi1.reserve(max_size*dim);
    nablaphi1.reserve(max_size*nabla_dim);
       
    vector <vector < adept::adouble> > vx(dim);
    vector <vector < adept::adouble> > vx_face(dim);
      
    for(int i=0;i<dim;i++){
      vx[i].reserve(max_size);
      vx_face[i].resize(9);
    }
   
    vector< vector< adept::adouble > > Soli(dim+1);
    vector< vector< int > > dofsVAR(dim+1); 
    for(int i=0;i<dim+1;i++){
      Soli[i].reserve(max_size);
      dofsVAR[i].reserve(max_size);
    }
    
    vector< vector< double > > Rhs(dim+1);
    vector< vector< adept::adouble > > aRhs(dim+1);
    for(int i=0;i<dim+1;i++){
      aRhs[i].reserve(max_size);
      Rhs[i].reserve(max_size);
    }     
    
    vector < int > dofsAll;
    dofsAll.reserve(max_size*(dim+1));
        
    vector < double > KKloc;
    KKloc.reserve(dim*max_size*(dim+1)*dim*max_size*(dim+1));
        
    vector < double > Jac;
    Jac.reserve(dim*max_size*(dim+1)*dim*max_size*(dim+1));
    
    // ------------------------------------------------------------------------
    // Physical parameters
    double rhof	 	= ml_prob.parameters.get<Fluid>("Fluid").get_density();             
    double IRe 		= ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number(); 
    double betans	= 1.;
       
    // gravity
    double _gravity[3]={0.,0.,0.};
     
    IRe =( 150*(counter+1) < 2000 )? 1./(150*(counter+1)):0.0005;
    cout<<"iteration="<<counter<<" Inverse Reynolds = "<<IRe<<endl;
    counter++;
    // -----------------------------------------------------------------
    // space discretization parameters
    unsigned SolType2 = ml_sol->GetSolutionType(ml_sol->GetIndex("U"));  
    unsigned end_ind2   = mymsh->GetEndIndex(SolType2);

    unsigned SolType1 = ml_sol->GetSolutionType(ml_sol->GetIndex("P"));  
    unsigned end_ind1   = mymsh->GetEndIndex(SolType1);

    // mesh and procs
    unsigned nel    = mymsh->GetElementNumber();
    unsigned igrid  = mymsh->GetGridNumber();
    unsigned iproc  = mymsh->processor_id();

    //----------------------------------------------------------------------------------
    //variable-name handling
    const char varname[4][3] = {"U","V","W","P"};
    vector <unsigned> indexVAR(dim+1);
    vector <unsigned> indVAR(dim+1);  
    vector <unsigned> SolType(dim+1);  
  
    for(unsigned ivar=0; ivar<dim; ivar++) {
      indVAR[ivar]=ml_sol->GetIndex(&varname[ivar][0]);
      SolType[ivar]=ml_sol->GetSolutionType(&varname[ivar][0]);
      indexVAR[ivar]=my_nnlin_impl_sys.GetSolPdeIndex(&varname[ivar][0]);  
    }
    indVAR[dim]=ml_sol->GetIndex(&varname[3][0]);
    SolType[dim]=ml_sol->GetSolutionType(&varname[3][0]);
    indexVAR[dim]=my_nnlin_impl_sys.GetSolPdeIndex(&varname[3][0]);  
    
    //----------------------------------------------------------------------------------
        
    start_time=clock();
    
    myKK->zero();
    
    // *** element loop ***
    for(int iel=mymsh->IS_Mts2Gmt_elem_offset[iproc]; iel < mymsh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

      unsigned kel        = mymsh->IS_Mts2Gmt_elem[iel]; 
      short unsigned kelt = myel->GetElementType(kel);
      unsigned nve        = myel->GetElementDofNumber(kel,end_ind2);
      unsigned nve1       = myel->GetElementDofNumber(kel,end_ind1);
     

      // *******************************************************************************************************
    
      //initialization of everything is in common fluid and solid
    
      //Rhs
      for(int i=0; i<dim; i++) {
	dofsVAR[i].resize(nve);
	Soli[indexVAR[i]].resize(nve);
	aRhs[indexVAR[i]].resize(nve);
	Rhs[indexVAR[i]].resize(nve);
      }
      dofsVAR[dim].resize(nve1);
      Soli[indexVAR[dim]].resize(nve1);
      aRhs[indexVAR[dim]].resize(nve1);
      Rhs[indexVAR[dim]].resize(nve1);
      
      dofsAll.resize(0);
      
      KKloc.resize((dim*nve+nve1)*(dim*nve+nve1));
      Jac.resize((dim*nve+nve1)*(dim*nve+nve1));
      
      // ----------------------------------------------------------------------------------------
      // coordinates, solutions, displacement, velocity dofs
        
      for(int i=0;i<dim;i++){
	vx[i].resize(nve);
      }
    
      for (unsigned i=0;i<nve;i++) {
	unsigned inode=myel->GetMeshDof(kel,i,SolType2);
	unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
	// flag to know if the node "inode" lays on the fluid-solid interface
	
	for(int j=0; j<dim; j++) {
	  // velocity dofs
	  Soli[indexVAR[j]][i] =  (*mysolution->_Sol[indVAR[j]])(inode_Metis);
	  aRhs[indexVAR[j]][i] = 0.;
	  dofsVAR[j][i] = myLinEqSolver->GetKKDof(indVAR[j],indexVAR[j],inode); 
	  //coordinates
	  vx[j][i]=  (*mymsh->_coordinate->_Sol[j])(inode_Metis); 
	  
	}
      }

      // pressure dofs
      for (unsigned i=0;i<nve1;i++) {
	unsigned inode=myel->GetMeshDof(kel,i,SolType1);
	unsigned inode_Metis =mymsh->GetMetisDof(inode,SolType[dim]);
	dofsVAR[dim][i]=myLinEqSolver->GetKKDof(indVAR[dim],indexVAR[dim],inode);
	Soli[indexVAR[dim]][i] = (*mysolution->_Sol[indVAR[dim]])(inode_Metis);
	aRhs[indexVAR[dim]][i] = 0.;
      }
      
      // build dof ccomposition             
      for(int idim=0;idim<dim;idim++){
	dofsAll.insert( dofsAll.end(), dofsVAR[idim].begin(), dofsVAR[idim].end() );
      }
      dofsAll.insert( dofsAll.end(), dofsVAR[dim].begin(), dofsVAR[dim].end() );
 
      if (igrid==gridn || !myel->GetRefinedElementIndex(kel) ) {  
	
	adeptStack.new_recording();
		
// 	// Boundary integral
// 	{
// 	  double tau=0.;
// 	  vector<adept::adouble> normal(dim,0);
// 	       
// 	  // loop on faces
// 	  for(unsigned jface=0; jface<myel->GetElementFaceNumber(kel); jface++) {
// 		
// 	    // look for boundary faces
// 	    if(myel->GetFaceElementIndex(kel,jface)<0) {
// 	      unsigned int face = -(mymsh->el->GetFaceElementIndex(kel,jface)+1);	      
// 	      if( !ml_sol->_SetBoundaryConditionFunction(0.,0.,0.,"U",tau,face,0.) && tau!=0.){
// 		unsigned nve = mymsh->el->GetElementFaceDofNumber(kel,jface,SolType2);
// 		const unsigned felt = mymsh->el->GetElementFaceType(kel, jface);  		  		  
// 		for(unsigned i=0; i<nve; i++) {
// 		  unsigned inode=mymsh->el->GetFaceVertexIndex(kel,jface,i)-1u;
// 		  unsigned inode_Metis=mymsh->GetMetisDof(inode,2);
// 		  unsigned int ilocal = mymsh->el->GetLocalFaceVertexIndex(kel, jface, i);
// 		  for(unsigned idim=0; idim<dim; idim++) {
// 		    vx_face[idim][i]=(*mymsh->_coordinate->_Sol[idim])(inode_Metis) + Soli[indexVAR[idim]][ilocal];
// 		  }
// 		}
// 		for(unsigned igs=0; igs < mymsh->_finiteElement[felt][SolType2]->GetGaussPointNumber(); igs++) {
// 		  mymsh->_finiteElement[felt][SolType2]->JacobianSur_AD(vx_face,igs,Weight,phi,gradphi,normal);
// 		  //phi1 =mymsh->_finiteElement[felt][SolType2]->GetPhi(igs);
// 		  // *** phi_i loop ***
// 		  for(unsigned i=0; i<nve; i++) {
// 		    adept::adouble value = - phi[i]*tau/rhof*Weight;
// 		    unsigned int ilocal = mymsh->el->GetLocalFaceVertexIndex(kel, jface, i);
// 		    
// 		    for(unsigned idim=0; idim<dim; idim++) {
// 		      if(1){
// 			aRhs[indexVAR[dim+idim]][ilocal]   += value*normal[idim];
// 		      }
// 		      else { //if interface node it goes to solid
// 			aRhs[indexVAR[idim]][ilocal]   += value*normal[idim];
// 		      }
// 		    }	    
// 		  }
// 		}
// 	      }
// 	    }
// 	  }    
// 	}
	  	  
	// *** Gauss point loop ***
	
	//adept::adouble supg_tau;
	for (unsigned ig=0;ig < mymsh->_finiteElement[kelt][SolType2]->GetGaussPointNumber(); ig++) {
	  // *** get Jacobian and test function and test function derivatives in the moving frame***
	  mymsh->_finiteElement[kelt][SolType2]->Jacobian_AD(vx,ig,Weight,phi,gradphi,nablaphi);
	  mymsh->_finiteElement[kelt][SolType1]->Jacobian_AD(vx,ig,Weight1,phi1,gradphi1,nablaphi1);  
	  
	 
	 // vector< double > V(dim,0.);
	  //unsigned ir = referenceElementPoint[kelt];
	
	  //  velocity: solution, gradient and laplace
	  for(int i=0; i<dim; i++){
	    SolVAR[i]=0.;
	    for(int j=0; j<dim; j++) {
	      GradSolVAR[i][j]=0.;
	      NablaSolVAR[i][j]=0.;
	    }
	    for (unsigned inode=0; inode<nve; inode++) {
	      adept::adouble soli = Soli[indexVAR[i]][inode];
	      //if(inode==ir) V[i]=soli.value();
	      SolVAR[i]+=phi[inode]*soli;
	      for(int j=0; j<dim; j++) {
		GradSolVAR[i][j]+=gradphi[inode*dim+j]*soli;
		NablaSolVAR[i][j]+=nablaphi[inode*nabla_dim+j]*soli;
	      }
	    }
	  } 
	  // pressure, solution and gradient 
	  SolVAR[dim]=0.;
	  for(int j=0; j<dim; j++) {
	      GradSolVAR[dim][j]=0.;
	  }
	  for (unsigned inode=0; inode<nve1; inode++) {
	    adept::adouble soli = Soli[indexVAR[dim]][inode];
	    SolVAR[dim]+=phi1[inode]*soli;
	    for(int j=0; j<dim; j++) {
	      GradSolVAR[dim][j]+=gradphi1[inode*dim+j]*soli;
	    }
	  } 
	  
	  //BEGIN TAU_SUPG EVALUATION ============
	  // ********************************* Tau_Supg ******************************************
	  // Computer Methods in Applied Mechanics and Engineering 95 (1992) 221-242 North-Holland
	  // *************************************************************************************	  
	  // velocity
	  vector < double > u(dim);
	  for(int ivar=0; ivar<dim; ivar++){
	    u[ivar]=SolVAR[ivar].value();
	  }
	  
	  // speed
	  double uL2Norm=0.;
	  for(int ivar=0;ivar<dim;ivar++){
	    uL2Norm += u[ivar]*u[ivar];
	  }
	  uL2Norm=sqrt(uL2Norm);
	  double tauSupg=0.;
	  if( uL2Norm/(2.*IRe) > 1.0e-10){
	    // velocity direction s = u/|u|
	    vector < double > s(dim);
	    for(int ivar=0;ivar<dim;ivar++)
	      s[ivar]=u[ivar]/uL2Norm;
	  
	    // element lenght h(s) = 2. ( \sum_i |s . gradphi_i | )^(-1)
	    double h=0;
	    for (unsigned i=0; i<nve; i++) {
	      double sDotGradphi=0.; 
	      for(int ivar=0; ivar<dim; ivar++)
		sDotGradphi += s[ivar]*gradphi[i*dim+ivar].value();
	      h += fabs(sDotGradphi);
	    }
	    h = 2./h;
	  
	    //tauSupg
	    double Reu   = (uL2Norm*h)/(2*IRe);
	    double zReu  = (Reu <= 3)? Reu/3.:1;
	    tauSupg = h / (2.*uL2Norm)*zReu;
	  }
	  //END TAU_SUPG EVALUATION ============
	  
	  //BEGIN FLUID ASSEMBLY ============
	  { 
	    vector < adept::adouble > ResSupg(dim,0.);
	    for(unsigned ivar=0; ivar<dim; ivar++) {
	      for(unsigned jvar=0; jvar<dim; jvar++) {
		ResSupg[ivar] += SolVAR[jvar]*GradSolVAR[ivar][jvar] - IRe*NablaSolVAR[ivar][jvar];
	      }
	      ResSupg[ivar] += GradSolVAR[dim][ivar];
	    }
	    
	    //BEGIN redidual momentum block  
	    for (unsigned i=0; i<nve; i++){
   
	      
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		adept::adouble Adv_rhs=0.;
		adept::adouble Lap_rhs=0.;
		adept::adouble supgPhi=0.;
		for(unsigned jvar=0; jvar<dim; jvar++) {
		  Lap_rhs += gradphi[i*dim+jvar]*GradSolVAR[ivar][jvar];
		  Adv_rhs += SolVAR[jvar]*GradSolVAR[ivar][jvar];
		  supgPhi += SolVAR[jvar]*gradphi[i*dim+jvar]*tauSupg; 
		}
		aRhs[indexVAR[ivar]][i]+= ( -IRe*Lap_rhs - Adv_rhs*phi[i] + SolVAR[dim]*gradphi[i*dim+ivar]
					    - ResSupg[ivar]*supgPhi )*Weight;	      
	      }
	    } 
	    //END redidual momentum block     
	    
	    //BEGIN continuity block 
	    {  	    
	      adept::adouble div_vel=0.;
	      for(int i=0; i<dim; i++) {
		div_vel +=GradSolVAR[i][i];
	      }
	      for (unsigned i=0; i<nve1; i++) {
		
		adept::adouble supgPhi=0.;
		for(int ivar=0;ivar<dim;ivar++){
		  supgPhi += SolVAR[ivar]*gradphi1[i*dim+ivar]*tauSupg; 
		}
		aRhs[indexVAR[dim]][i] += -( (phi1[i]+0*supgPhi) * (-div_vel) )*Weight;
	      }
	    }
	    //END continuity block ===========================
	  }   
	  //END FLUID ASSEMBLY ============
	}
      }
	
      //BEGIN local to global assembly 	
      //copy adouble aRhs into double Rhs
      for (unsigned i=0;i<dim;i++) {
	for(int j=0; j<nve; j++) {
	  Rhs[indexVAR[i]][j] = aRhs[indexVAR[i]][j].value();
	}
      }
      for (unsigned j=0;j<nve1;j++) {
	Rhs[indexVAR[dim]][j] = aRhs[indexVAR[dim]][j].value();
      }	
      for(int i=0; i<dim+1; i++) {
	myRES->add_vector_blocked(Rhs[indexVAR[i]],dofsVAR[i]);
      }     

      //Store equations
      for(int i=0; i<dim; i++) {  
	adeptStack.dependent(&aRhs[indexVAR[i]][0], nve);
	adeptStack.independent(&Soli[indexVAR[i]][0], nve); 
      }
      adeptStack.dependent(&aRhs[indexVAR[dim]][0], nve1);
      adeptStack.independent(&Soli[indexVAR[dim]][0], nve1);   
      adeptStack.jacobian(&Jac[0]);	
      unsigned nveAll=(dim*nve+nve1);
      for (int inode=0;inode<nveAll;inode++){
	for (int jnode=0;jnode<nveAll;jnode++){
	   KKloc[inode*nveAll+jnode]=-Jac[jnode*nveAll+inode];
	}
      }
      myKK->add_matrix_blocked(KKloc,dofsAll,dofsAll);
      adeptStack.clear_independents();
      adeptStack.clear_dependents();
       
      //END local to global assembly
   
    } //end list of elements loop
    
    myKK->close();
    myRES->close();
    
    // *************************************
    end_time=clock();
    AssemblyTime+=(end_time-start_time);
    // ***************** END ASSEMBLY RESIDUAL + MATRIX *******************

  }  
  
  

void AssembleMatrixResNS_old(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix){
     
  //pointers 
  Solution*	 mysolution  	             = ml_prob._ml_sol->GetSolutionLevel(level);
  NonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<NonLinearImplicitSystem>("Navier-Stokes");
  LinearEquationSolver*  mylsyspde	     = my_nnlin_impl_sys._LinSolver[level];   
  const char* pdename                        = my_nnlin_impl_sys.name().c_str();
  
  MultiLevelSolution* ml_sol=ml_prob._ml_sol;
  
  
  mesh*		 mymsh    	= ml_prob._ml_msh->GetLevel(level);
  elem*		 myel		= mymsh->el;
  SparseMatrix*	 myKK		= mylsyspde->_KK;
  NumericVector* myRES 		= mylsyspde->_RES;
    
  //data
  const unsigned dim = mymsh->GetDimension();
  const unsigned nabla_dim = 3*(dim-1);
  unsigned nel= mymsh->GetElementNumber();
  unsigned igrid= mymsh->GetGridNumber();
  unsigned iproc = mymsh->processor_id();
  double ILambda= 0; 
  double IRe = ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
  bool penalty = true; 
  const bool symm_mat = false;
  const bool NavierStokes = true; 
  unsigned nwtn_alg = 2; 
  bool newton = (nwtn_alg==0) ? 0:1; 
  
  
  if(counter==0) IRe=0.1;
  counter++;
  
  // solution and coordinate variables
  const char Solname[4][2] = {"U","V","W","P"};
  vector < unsigned > SolPdeIndex(dim+1);
  vector < unsigned > SolIndex(dim+1);  
  
  //const char coordinate_name[3][2] = {"X","Y","Z"};
  //vector < unsigned > coordinate_Index(dim);
  vector< vector < double> > coordinates(dim);
  
  for(unsigned ivar=0; ivar<dim; ivar++) {
    SolPdeIndex[ivar]=my_nnlin_impl_sys.GetSolPdeIndex(&Solname[ivar][0]);
    SolIndex[ivar]=ml_sol->GetIndex(&Solname[ivar][0]);
    //coordinate_Index[ivar]=ivar;//ml_prob.GetIndex(&coordinate_name[ivar][0]);
  }
  SolPdeIndex[dim]=my_nnlin_impl_sys.GetSolPdeIndex(&Solname[3][0]);
  SolIndex[dim]=ml_sol->GetIndex(&Solname[3][0]);       
  //solution order
  unsigned order_ind_vel = ml_sol->GetSolutionType(SolIndex[0]);
  unsigned end_ind2   = mymsh->GetEndIndex(order_ind_vel);
  unsigned order_ind_p = ml_sol->GetSolutionType(SolIndex[dim]);
  unsigned end_ind1   = mymsh->GetEndIndex(order_ind_p);
  
//   double alpha = 0.;
//   if(order_ind_p == order_ind_vel && order_ind_vel == 0) // if pressure and velocity are both linear, we need stabilization 
//   {
//     alpha = 0.005; 
//   }
  
  // declare 
  vector < int > metis_node2; 
  vector < int > node1;
  vector< vector< int > > KK_dof(dim+1); 
  vector <double> phi2;
  vector <double> gradphi2;
  vector <double> nablaphi2;
  double Weight2;
  
  vector <double> phi1;
  vector <double> gradphi1;
  vector <double> nablaphi1;
  double Weight1;
  
  double normal[3];
  vector< vector< double > > F(dim+1);
  vector< vector< vector< double > > > B(dim+1); 
  
  // reserve
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
  metis_node2.reserve(max_size);
  node1.reserve( static_cast< unsigned > (ceil(pow(2,dim))));
  for(int i=0;i<dim;i++) {
    coordinates[i].reserve(max_size);
  }
  phi2.reserve(max_size);
  gradphi2.reserve(max_size*dim);	
  nablaphi2.reserve(max_size*nabla_dim);
  
  phi1.reserve(max_size);
  gradphi1.reserve(max_size*dim);	
  nablaphi1.reserve(max_size*nabla_dim);
  
  for(int i=0;i<dim;i++) {
    KK_dof[i].reserve(max_size);
  }
   
  for(int i=0;i<dim+1;i++) F[i].reserve(max_size);
    
  if(assembe_matrix){
    for(int i=0;i<dim+1;i++){
      B[i].resize(dim+1);
      for(int j=0;j<dim+1;j++){
	B[i][j].reserve(max_size*max_size);
      }
    }
  }
    
  vector < double > SolVAR(dim+1);
  vector < vector < double > > gradSolVAR(dim);
  vector < vector < double > > NablaSolVAR(dim);
  for(int i=0;i<dim;i++) {
    gradSolVAR[i].resize(dim);  
    NablaSolVAR[i].resize(dim); 
  }
  
  // Set to zeto all the entries of the matrix
  if(assembe_matrix) myKK->zero();
  
  // *** element loop ***
 
  for (int iel=mymsh->IS_Mts2Gmt_elem_offset[iproc]; iel < mymsh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

    unsigned kel = mymsh->IS_Mts2Gmt_elem[iel];
    short unsigned kelt=myel->GetElementType(kel);
    unsigned nve2=myel->GetElementDofNumber(kel,end_ind2);
    unsigned nve1=myel->GetElementDofNumber(kel,end_ind1);
    
    //set to zero all the entries of the FE matrices
    metis_node2.resize(nve2);
    node1.resize(nve1);
    
    for(int ivar=0; ivar<dim; ivar++) {
      coordinates[ivar].resize(nve2);
      KK_dof[ivar].resize(nve2);
      
      F[SolPdeIndex[ivar]].resize(nve2);
      memset(&F[SolPdeIndex[ivar]][0],0,nve2*sizeof(double));
      
      if(assembe_matrix){
	B[SolPdeIndex[ivar]][SolPdeIndex[ivar]].resize(nve2*nve2);
	B[SolPdeIndex[ivar]][SolPdeIndex[dim]].resize(nve2*nve1);
	B[SolPdeIndex[dim]][SolPdeIndex[ivar]].resize(nve1*nve2);
	memset(&B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][0],0,nve2*nve2*sizeof(double));
	memset(&B[SolPdeIndex[ivar]][SolPdeIndex[dim]][0],0,nve2*nve1*sizeof(double));
	memset(&B[SolPdeIndex[dim]][SolPdeIndex[ivar]][0],0,nve1*nve2*sizeof(double));
      }
    }
    KK_dof[dim].resize(nve1);
    F[SolPdeIndex[dim]].resize(nve1);
    memset(&F[SolPdeIndex[dim]][0],0,nve1*sizeof(double));
      
      
    if(assembe_matrix*nwtn_alg==2){
      for(int ivar=0; ivar<dim; ivar++) {
	for(int jvar=1; jvar<dim; jvar++) {
	  B[SolPdeIndex[ivar]][SolPdeIndex[(ivar+jvar)%dim]].resize(nve2*nve2);
	  memset(&B[SolPdeIndex[ivar]][SolPdeIndex[(ivar+jvar)%dim]][0],0,nve2*nve2*sizeof(double));
	}
      }
    }
  
    if(assembe_matrix*penalty){
      B[SolPdeIndex[dim]][SolPdeIndex[dim]].resize(nve1*nve1,0.);
      memset(&B[SolPdeIndex[dim]][SolPdeIndex[dim]][0],0,nve1*nve1*sizeof(double));
    }
    
    for( unsigned i=0;i<nve2;i++){
      unsigned inode=myel->GetMeshDof(kel,i,order_ind_vel);
      //unsigned inode=myel->GetElementVertexIndex(kel,i)-1u;
      unsigned inode_coord_metis=mymsh->GetMetisDof(inode,2);
      metis_node2[i]=mymsh->GetMetisDof(inode,order_ind_vel);
      for(unsigned ivar=0; ivar<dim; ivar++) {
	coordinates[ivar][i]=(*mymsh->_coordinate->_Sol[ivar])(inode_coord_metis);
	KK_dof[ivar][i]=mylsyspde->GetKKDof(SolIndex[ivar],SolPdeIndex[ivar],inode);
      }
    }
    
    //double hk = sqrt( (coordinates[0][2] - coordinates[0][0])*(coordinates[0][2] - coordinates[0][0]) + 
     // (coordinates[1][2] - coordinates[1][0])*(coordinates[1][2] - coordinates[1][0]) );
    
    for(unsigned i=0;i<nve1;i++) {
      unsigned inode=myel->GetMeshDof(kel,i,order_ind_p);
      //unsigned inode=(order_ind_p<dim)?(myel->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
      node1[i]=inode;
      KK_dof[dim][i]=mylsyspde->GetKKDof(SolIndex[dim],SolPdeIndex[dim],inode);
    }
   
    if(igrid==gridn || !myel->GetRefinedElementIndex(kel)) {
      vector< double > V(dim,0.);
      unsigned ir = referenceElementPoint[kelt];
      
      // *** Gauss poit loop ***
      for(unsigned ig=0;ig < ml_prob._ml_msh->_finiteElement[kelt][order_ind_vel]->GetGaussPointNumber(); ig++) {
	// *** get Jacobian and test function and test function derivatives ***
	ml_prob._ml_msh->_finiteElement[kelt][order_ind_vel]->Jacobian(coordinates,ig,Weight2,phi2,gradphi2,nablaphi2);
	ml_prob._ml_msh->_finiteElement[kelt][order_ind_p]->Jacobian(coordinates,ig,Weight1,phi1,gradphi1,nablaphi1);

	double GradSolP[3] = {0.,0.,0.};
	//velocity variable
	for(unsigned ivar=0; ivar<dim; ivar++) {
	  SolVAR[ivar]=0;
	  for(unsigned jvar=0; jvar<dim; jvar++){ 
	    gradSolVAR[ivar][jvar]=0; 
	    NablaSolVAR[ivar][jvar]=0.; 
	  }
	  unsigned SolIndex=ml_sol->GetIndex(&Solname[ivar][0]);
	  unsigned SolType=ml_sol->GetSolutionType(&Solname[ivar][0]);
	  for(unsigned i=0; i<nve2; i++) {
	    double soli = (*mysolution->_Sol[SolIndex])(metis_node2[i]);
	    
	    if(i == ir) V[ivar]=soli;
	    
	    SolVAR[ivar]+=phi2[i]*soli;
	    for(unsigned jvar=0; jvar<dim; jvar++){
	      gradSolVAR[ivar][jvar]  += gradphi2[i*dim+jvar]*soli; 
	      NablaSolVAR[ivar][jvar] += nablaphi2[i*nabla_dim+jvar]*soli;
	    }
	  }
	}
	//pressure variable
	SolVAR[dim]=0;
	unsigned SolIndex=ml_sol->GetIndex(&Solname[3][0]);
	unsigned SolType=ml_sol->GetSolutionType(&Solname[3][0]);
	for(unsigned i=0; i<nve1; i++){
	  unsigned sol_dof = mymsh->GetMetisDof(node1[i],SolType);
	  double soli = (*mysolution->_Sol[SolIndex])(sol_dof);
	  SolVAR[dim]+=phi1[i]*soli;
	  for(unsigned ivar=0; ivar<dim; ivar++){
	    GradSolP[ivar] += gradphi1[i*dim+ivar]*soli;
	  }
	}
	
	// Supg stabilization tau evaluation
	
	//vector< double > V(dim,0.);
// 	for(int ivar=0;ivar<dim;ivar++) 
// 	  V[ivar] = SolVAR[ivar];
      	double nu=IRe;
	double barNu=0.;
	double vL2Norm2=0.;
	for(int ivar=0;ivar<dim;ivar++){
	  vL2Norm2 += V[ivar]*V[ivar];
	  unsigned ip = referenceElementDirection[kelt][ivar][1];
	  unsigned im = referenceElementDirection[kelt][ivar][0];
	  double VxiHxi=0.;
	  for(int j=0;j<dim;j++){
	    VxiHxi += (coordinates[j][ip]-coordinates[j][im]) * V[j];
	  }	
	  double PeXi=VxiHxi/(2.*nu);		
	  double barXi = ( fabs( PeXi ) < 1.0e-10) ? 0. : 1./tanh(PeXi)-1./PeXi;
	  barNu += barXi * VxiHxi /2.;
	}
	double supgTau = ( vL2Norm2 > 1.0e-15 ) ? 5*barNu/vL2Norm2 : 0.;
	// End Stabilization stabilization tau evaluation
	//supgTau=0.;
	
	// *** phi_i loop ***
	for(unsigned i=0; i<nve2; i++){
	
	  //BEGIN RESIDUALS A block ===========================
	  double supgPhi;
	  vector <double> Adv_rhs(dim,0.); 
	  for(unsigned ivar=0; ivar<dim; ivar++) {
	    Adv_rhs[ivar]=0.;
	    double Lap_rhs=0.;
	    //double resRhs=0.;
	    supgPhi=0.;
	    for(unsigned jvar=0; jvar<dim; jvar++) {
	      Lap_rhs += gradphi2[i*dim+jvar]*gradSolVAR[ivar][jvar];
	      Adv_rhs[ivar] += SolVAR[jvar]*gradSolVAR[ivar][jvar];
	      //resRhs  += SolVAR[jvar]*gradSolVAR[ivar][jvar]-0*IRe*NablaSolVAR[ivar][jvar];
	      supgPhi += (SolVAR[jvar]*gradphi2[i*dim+jvar])* supgTau; 
	    }
	    //resRhs += GradSolP[ivar];
	    F[SolPdeIndex[ivar]][i]+= ( -IRe*Lap_rhs-NavierStokes*Adv_rhs[ivar]*(phi2[i]+supgPhi)
					+SolVAR[dim]*gradphi2[i*dim+ivar])*Weight2;
	  }
	  //END RESIDUALS A block ===========================
	  
	  if(assembe_matrix){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve2; j++) {
	      double Lap=0;
	      
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		// Laplacian
		Lap  += (gradphi2[i*dim+ivar]*gradphi2[j*dim+ivar])*Weight2;
	      }
	      for(unsigned ivar=0; ivar<dim; ivar++) {    
		B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][i*nve2+j] += IRe*Lap;
		  		  		  
		  for(unsigned jvar=0; jvar<dim; jvar++) {
		    B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][i*nve2+j] += (phi2[i]+supgPhi)*SolVAR[jvar]*gradphi2[j*dim+jvar]*Weight2;
		    B[SolPdeIndex[ivar]][SolPdeIndex[jvar]][i*nve2+j] += (phi2[i]+supgPhi)*phi2[j]*gradSolVAR[ivar][jvar]*Weight2;
		    B[SolPdeIndex[ivar]][SolPdeIndex[jvar]][i*nve2+j] +=  phi2[j]*gradphi2[i*dim+jvar]*supgTau*Adv_rhs[ivar]*Weight2;
		  }
		  
		
	      }
  	    } //end phij loop
	    
	    // *** phi1_j loop ***
	    for(unsigned j=0; j<nve1; j++){
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndex[ivar]][SolPdeIndex[dim]][i*nve1+j] -= gradphi2[i*dim+ivar]*phi1[j]*Weight2;
	      }
	    } //end phi1_j loop
	  } // endif assembe_matrix
	} //end phii loop
  

	// *** phi1_i loop ***
	for(unsigned i=0; i<nve1; i++){
	  //BEGIN RESIDUALS B block ===========================
	  double div = 0;
	  for(unsigned ivar=0; ivar<dim; ivar++) {
	    div += gradSolVAR[ivar][ivar];
	  }
	  F[SolPdeIndex[dim]][i]+= (phi1[i]*div)*Weight2;// + /*penalty*ILambda*phi1[i]*SolVAR[dim]*/ 
	                             //+ 0*hk*hk*(1./0.001)*alpha*(GradSolP[0]*gradphi2[i*dim + 0] + GradSolP[1]*gradphi2[i*dim + 1]) )*Weight2;
				     
           
	  //END RESIDUALS  B block ===========================
	  
	  if(assembe_matrix){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve2; j++) {
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndex[dim]][SolPdeIndex[ivar]][i*nve2+j]-= phi1[i]*gradphi2[j*dim+ivar]*Weight2;
	      }
	    }  //end phij loop
	  } // endif assembe_matrix
	}  //end phi1_i loop
	
	if(assembe_matrix * penalty){  //block nve1 nve1
	  // *** phi_i loop ***
	  for(unsigned i=0; i<nve1; i++){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve1; j++){
	      //B[SolPdeIndex[dim]][SolPdeIndex[dim]][i*nve1+j]-= ILambda*phi1[i]*phi1[j]*Weight2;
	      for(unsigned ivar=0; ivar<dim; ivar++) {
	        B[SolPdeIndex[dim]][SolPdeIndex[dim]][i*nve1+j] -= 0.;//*hk*hk*(1./0.001)*alpha*(gradphi2[i*dim + ivar]*gradphi2[j*dim + ivar])*Weight2; 
	      }
	    }
	  }
	}   //end if penalty
      }  // end gauss point loop
      
      //--------------------------------------------------------------------------------------------------------
      // Boundary Integral --> to be added
      //number of faces for each type of element
//       if (igrid==gridn || !myel->GetRefinedElementIndex(kel) ) {
//      
// 	unsigned nfaces = myel->GetElementFaceNumber(kel);
// 
// 	// loop on faces
// 	for(unsigned jface=0;jface<nfaces;jface++){ 
// 	  
// 	  // look for boundary faces
// 	  if(myel->GetFaceElementIndex(kel,jface)<0){
// 	    for(unsigned ivar=0; ivar<dim; ivar++) {
// 	      ml_prob.ComputeBdIntegral(pdename, &Solname[ivar][0], kel, jface, level, ivar);
// 	    }
// 	  }
// 	}	
//       }
      //--------------------------------------------------------------------------------------------------------
    } // endif single element not refined or fine grid loop
    //--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the Global Matrix/Vector
    for(unsigned ivar=0; ivar<dim; ivar++) {
      myRES->add_vector_blocked(F[SolPdeIndex[ivar]],KK_dof[ivar]);
      if(assembe_matrix){
	myKK->add_matrix_blocked(B[SolPdeIndex[ivar]][SolPdeIndex[ivar]],KK_dof[ivar],KK_dof[ivar]);  
	myKK->add_matrix_blocked(B[SolPdeIndex[ivar]][SolPdeIndex[dim]],KK_dof[ivar],KK_dof[dim]);
	myKK->add_matrix_blocked(B[SolPdeIndex[dim]][SolPdeIndex[ivar]],KK_dof[dim],KK_dof[ivar]);
	if(nwtn_alg==2){
	  for(unsigned jvar=1; jvar<dim; jvar++) {
	    myKK->add_matrix_blocked(B[SolPdeIndex[ivar]][SolPdeIndex[(ivar+jvar)%dim]],KK_dof[ivar],KK_dof[(ivar+jvar)%dim]);  
	  }
	}
      }
    }
    //Penalty
    if(assembe_matrix*penalty) myKK->add_matrix_blocked(B[SolPdeIndex[dim]][SolPdeIndex[dim]],KK_dof[dim],KK_dof[dim]);
    myRES->add_vector_blocked(F[SolPdeIndex[dim]],KK_dof[dim]);
    //--------------------------------------------------------------------------------------------------------  
  } //end list of elements loop for each subdomain
  
  
  if(assembe_matrix) myKK->close();
  myRES->close();
  // ***************** END ASSEMBLY *******************
}

  
  
  
  
  
  
  /*
  
  
  
  
  
  
  static adept::Stack s; 
  
  //pointers 
  Solution*	 mysolution  	             = ml_prob._ml_sol->GetSolutionLevel(level);
  NonLinearImplicitSystem& my_nnlin_impl_sys = ml_prob.get_system<NonLinearImplicitSystem>("Navier-Stokes");
  LinearEquationSolver*  mylsyspde	     = my_nnlin_impl_sys._LinSolver[level];   
  const char* pdename                        = my_nnlin_impl_sys.name().c_str();
  
  MultiLevelSolution* ml_sol=ml_prob._ml_sol;
  
  
  mesh*		 mymsh    	= ml_prob._ml_msh->GetLevel(level);
  elem*		 myel		= mymsh->el;
  SparseMatrix*	 myKK		= mylsyspde->_KK;
  NumericVector* myRES 		= mylsyspde->_RES;
    
  //data
  const unsigned dim = mymsh->GetDimension();
  const unsigned nabla_dim = 3*(dim-1);
  unsigned nel= mymsh->GetElementNumber();
  unsigned igrid= mymsh->GetGridNumber();
  unsigned iproc = mymsh->processor_id();
  double ILambda= 0; 
  double IRe = ml_prob.parameters.get<Fluid>("Fluid").get_IReynolds_number();
  bool penalty = true; 
  const bool symm_mat = false;
  const bool NavierStokes = true; 
  unsigned nwtn_alg = 2; 
  bool newton = (nwtn_alg==0) ? 0:1; 
  
  
  if(counter==0) IRe=0.1;
  counter++;
  
  // solution and coordinate variables
  const char Solname[4][2] = {"U","V","W","P"};
  vector < unsigned > SolPdeIndex(dim+1);
  vector < unsigned > SolIndex(dim+1);  
  
  //const char coordinate_name[3][2] = {"X","Y","Z"};
  //vector < unsigned > coordinate_Index(dim);
  vector< vector < double> > coordinates(dim);
  vector< vector < adept::adouble> > acoordinates(dim);
  
  for(unsigned ivar=0; ivar<dim; ivar++) {
    SolPdeIndex[ivar]=my_nnlin_impl_sys.GetSolPdeIndex(&Solname[ivar][0]);
    SolIndex[ivar]=ml_sol->GetIndex(&Solname[ivar][0]);
    //coordinate_Index[ivar]=ivar;//ml_prob.GetIndex(&coordinate_name[ivar][0]);
  }
  SolPdeIndex[dim]=my_nnlin_impl_sys.GetSolPdeIndex(&Solname[3][0]);
  SolIndex[dim]=ml_sol->GetIndex(&Solname[3][0]);       
  //solution order
  unsigned order_ind_vel = ml_sol->GetSolutionType(SolIndex[0]);
  unsigned end_ind2   = mymsh->GetEndIndex(order_ind_vel);
  unsigned order_ind_p = ml_sol->GetSolutionType(SolIndex[dim]);
  unsigned end_ind1   = mymsh->GetEndIndex(order_ind_p);
  
  // declare 
  vector < int > metis_node2; 
  vector < int > node1;
  vector< vector< int > > KK_dof(dim+1); 
  vector <double> phi2;
  vector <double> gradphi2;
  vector <double> nablaphi2;
  double Weight2;
  
  vector <adept::adouble> agradphi2;
  vector <adept::adouble> anablaphi2;
  adept::adouble aWeight2;
  
  
  
  vector <double> phi1;
  vector <double> gradphi1;
  vector <double> nablaphi1;
  double Weight1;
  
  vector <adept::adouble> agradphi1;
  vector <adept::adouble> anablaphi1;
  adept::adouble aWeight1;
  
  
  double normal[3];
  vector< vector< double > > F(dim+1);
  vector< vector< vector< double > > > B(dim+1); 
  
  // reserve
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
  metis_node2.reserve(max_size);
  node1.reserve( static_cast< unsigned > (ceil(pow(2,dim))));
  for(int i=0;i<dim;i++) {
    coordinates[i].reserve(max_size);
    acoordinates[i].reserve(max_size);
  }
  phi2.reserve(max_size);
  gradphi2.reserve(max_size*dim);	
  nablaphi2.reserve(max_size*nabla_dim);
  
  agradphi2.reserve(max_size*dim);	
  anablaphi2.reserve(max_size*nabla_dim);
  
  phi1.reserve(max_size);
  
  gradphi1.reserve(max_size*dim);	
  nablaphi1.reserve(max_size*nabla_dim);
  
  agradphi1.reserve(max_size*dim);	
  anablaphi1.reserve(max_size*nabla_dim);
  
  
  for(int i=0;i<dim;i++) {
    KK_dof[i].reserve(max_size);
  }
   
  for(int i=0;i<dim+1;i++) F[i].reserve(max_size);
    
  if(assembe_matrix){
    for(int i=0;i<dim+1;i++){
      B[i].resize(dim+1);
      for(int j=0;j<dim+1;j++){
	B[i][j].reserve(max_size*max_size);
      }
    }
  }
    
  // local objects
  vector<adept::adouble> aSolVAR(dim+1);
  vector<vector<adept::adouble> > aGradSolVAR(dim+1);
  vector<vector<adept::adouble> > aNablaSolVAR(dim+1);
      
  for(int i=0;i<dim+1;i++){
    aGradSolVAR[i].resize(dim);
    aNablaSolVAR[i].resize(dim);  
  }  
    
    
  //******************************************** 
  vector < double > SolVAR(dim+1);
  vector < vector < double > > gradSolVAR(dim);
  vector < vector < double > > NablaSolVAR(dim);
  for(int i=0;i<dim;i++) {
    gradSolVAR[i].resize(dim);  
    NablaSolVAR[i].resize(dim); 
  }
  //*********************************************
    
  vector< vector< adept::adouble > > Soli(dim+1);
  vector< vector< int > > dofsVAR(dim+1); 
  for(int i=0;i<dim+1;i++){
    Soli[i].reserve(max_size);
    dofsVAR[i].reserve(max_size);
  }
    
  vector< vector< double > > Rhs(dim+1);
  vector< vector< adept::adouble > > aRhs(dim+1);
  for(int i=0;i<dim+1;i++){
    aRhs[i].reserve(max_size);
    Rhs[i].reserve(max_size);
  }     
    
    
  vector < int > dofsAll;
  dofsAll.reserve(max_size*(dim+1));
        
  vector < double > KKloc;
  KKloc.reserve(dim*max_size*(dim+1)*dim*max_size*(dim+1));
        
  vector < double > Jac;
  Jac.reserve(dim*max_size*(dim+1)*dim*max_size*(dim+1));
    
  
  // **********************************************************
  
  
  
  
  
  
  
  
  
  
  
  // Set to zeto all the entries of the matrix
  if(assembe_matrix) myKK->zero();
  
  // *** element loop ***
 
  for (int iel=mymsh->IS_Mts2Gmt_elem_offset[iproc]; iel < mymsh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

    unsigned kel = mymsh->IS_Mts2Gmt_elem[iel];
    short unsigned kelt=myel->GetElementType(kel);
    unsigned nve2=myel->GetElementDofNumber(kel,end_ind2);
    unsigned nve1=myel->GetElementDofNumber(kel,end_ind1);
    
    //set to zero all the entries of the FE matrices
    metis_node2.resize(nve2);
    node1.resize(nve1);
    
    for(int ivar=0; ivar<dim; ivar++) {
      coordinates[ivar].resize(nve2);
      KK_dof[ivar].resize(nve2);
      
      F[SolPdeIndex[ivar]].resize(nve2);
      memset(&F[SolPdeIndex[ivar]][0],0,nve2*sizeof(double));
      
      if(assembe_matrix){
	B[SolPdeIndex[ivar]][SolPdeIndex[ivar]].resize(nve2*nve2);
	B[SolPdeIndex[ivar]][SolPdeIndex[dim]].resize(nve2*nve1);
	B[SolPdeIndex[dim]][SolPdeIndex[ivar]].resize(nve1*nve2);
	memset(&B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][0],0,nve2*nve2*sizeof(double));
	memset(&B[SolPdeIndex[ivar]][SolPdeIndex[dim]][0],0,nve2*nve1*sizeof(double));
	memset(&B[SolPdeIndex[dim]][SolPdeIndex[ivar]][0],0,nve1*nve2*sizeof(double));
      }
    }
    KK_dof[dim].resize(nve1);
    F[SolPdeIndex[dim]].resize(nve1);
    memset(&F[SolPdeIndex[dim]][0],0,nve1*sizeof(double));
      
      
    if(assembe_matrix*nwtn_alg==2){
      for(int ivar=0; ivar<dim; ivar++) {
	for(int jvar=1; jvar<dim; jvar++) {
	  B[SolPdeIndex[ivar]][SolPdeIndex[(ivar+jvar)%dim]].resize(nve2*nve2);
	  memset(&B[SolPdeIndex[ivar]][SolPdeIndex[(ivar+jvar)%dim]][0],0,nve2*nve2*sizeof(double));
	}
      }
    }
  
    if(assembe_matrix*penalty){
      B[SolPdeIndex[dim]][SolPdeIndex[dim]].resize(nve1*nve1,0.);
      memset(&B[SolPdeIndex[dim]][SolPdeIndex[dim]][0],0,nve1*nve1*sizeof(double));
    }
    
    for( unsigned i=0;i<nve2;i++){
      unsigned inode=myel->GetMeshDof(kel,i,order_ind_vel);
      //unsigned inode=myel->GetElementVertexIndex(kel,i)-1u;
      unsigned inode_coord_metis=mymsh->GetMetisDof(inode,2);
      metis_node2[i]=mymsh->GetMetisDof(inode,order_ind_vel);
      for(unsigned ivar=0; ivar<dim; ivar++) {
	coordinates[ivar][i]=(*mymsh->_coordinate->_Sol[ivar])(inode_coord_metis);
	KK_dof[ivar][i]=mylsyspde->GetKKDof(SolIndex[ivar],SolPdeIndex[ivar],inode);
      }
    }
    
    //double hk = sqrt( (coordinates[0][2] - coordinates[0][0])*(coordinates[0][2] - coordinates[0][0]) + 
     // (coordinates[1][2] - coordinates[1][0])*(coordinates[1][2] - coordinates[1][0]) );
    
    for(unsigned i=0;i<nve1;i++) {
      unsigned inode=myel->GetMeshDof(kel,i,order_ind_p);
      //unsigned inode=(order_ind_p<dim)?(myel->GetElementVertexIndex(kel,i)-1u):(kel+i*nel);
      node1[i]=inode;
      KK_dof[dim][i]=mylsyspde->GetKKDof(SolIndex[dim],SolPdeIndex[dim],inode);
    }
   
    if(igrid==gridn || !myel->GetRefinedElementIndex(kel)) {
      vector< double > V(dim,0.);
      unsigned ir = referenceElementPoint[kelt];
      
      // *** Gauss poit loop ***
      for(unsigned ig=0;ig < ml_prob._ml_msh->_finiteElement[kelt][order_ind_vel]->GetGaussPointNumber(); ig++) {
	// *** get Jacobian and test function and test function derivatives ***
	ml_prob._ml_msh->_finiteElement[kelt][order_ind_vel]->Jacobian(coordinates,ig,Weight2,phi2,gradphi2,nablaphi2);
	ml_prob._ml_msh->_finiteElement[kelt][order_ind_p]->Jacobian(coordinates,ig,Weight1,phi1,gradphi1,nablaphi1);

	double GradSolP[3] = {0.,0.,0.};
	//velocity variable
	for(unsigned ivar=0; ivar<dim; ivar++) {
	  SolVAR[ivar]=0;
	  for(unsigned jvar=0; jvar<dim; jvar++){ 
	    gradSolVAR[ivar][jvar]=0; 
	    NablaSolVAR[ivar][jvar]=0.; 
	  }
	  unsigned SolIndex=ml_sol->GetIndex(&Solname[ivar][0]);
	  unsigned SolType=ml_sol->GetSolutionType(&Solname[ivar][0]);
	  for(unsigned i=0; i<nve2; i++) {
	    double soli = (*mysolution->_Sol[SolIndex])(metis_node2[i]);
	    
	    if(i == ir) V[ivar]=soli;
	    
	    SolVAR[ivar]+=phi2[i]*soli;
	    for(unsigned jvar=0; jvar<dim; jvar++){
	      gradSolVAR[ivar][jvar]  += gradphi2[i*dim+jvar]*soli; 
	      NablaSolVAR[ivar][jvar] += nablaphi2[i*nabla_dim+jvar]*soli;
	    }
	  }
	}
	//pressure variable
	SolVAR[dim]=0;
	unsigned SolIndex=ml_sol->GetIndex(&Solname[3][0]);
	unsigned SolType=ml_sol->GetSolutionType(&Solname[3][0]);
	for(unsigned i=0; i<nve1; i++){
	  unsigned sol_dof = mymsh->GetMetisDof(node1[i],SolType);
	  double soli = (*mysolution->_Sol[SolIndex])(sol_dof);
	  SolVAR[dim]+=phi1[i]*soli;
	  for(unsigned ivar=0; ivar<dim; ivar++){
	    GradSolP[ivar] += gradphi1[i*dim+ivar]*soli;
	  }
	}
	
	// Supg stabilization tau evaluation
	
	//vector< double > V(dim,0.);
// 	for(int ivar=0;ivar<dim;ivar++) 
// 	  V[ivar] = SolVAR[ivar];
      	double nu=IRe;
	double barNu=0.;
	double vL2Norm2=0.;
	for(int ivar=0;ivar<dim;ivar++){
	  vL2Norm2 += V[ivar]*V[ivar];
	  unsigned ip = referenceElementDirection[kelt][ivar][1];
	  unsigned im = referenceElementDirection[kelt][ivar][0];
	  double VxiHxi=0.;
	  for(int j=0;j<dim;j++){
	    VxiHxi += (coordinates[j][ip]-coordinates[j][im]) * V[j];
	  }	
	  double PeXi=VxiHxi/(2.*nu);		
	  double barXi = ( fabs( PeXi ) < 1.0e-10) ? 0. : 1./tanh(PeXi)-1./PeXi;
	  barNu += barXi * VxiHxi /2.;
	}
	double supgTau = ( vL2Norm2 > 1.0e-15 ) ? 5*barNu/vL2Norm2 : 0.;
	// End Stabilization stabilization tau evaluation
	
	
	// *** phi_i loop ***
	for(unsigned i=0; i<nve2; i++){
	
	  //BEGIN RESIDUALS A block ===========================
	  double supgPhi;
	  vector <double> Adv_rhs(dim,0.); 
	  for(unsigned ivar=0; ivar<dim; ivar++) {
	    Adv_rhs[ivar]=0.;
	    double Lap_rhs=0.;
	    //double resRhs=0.;
	    supgPhi=0.;
	    for(unsigned jvar=0; jvar<dim; jvar++) {
	      Lap_rhs += gradphi2[i*dim+jvar]*gradSolVAR[ivar][jvar];
	      Adv_rhs[ivar] += SolVAR[jvar]*gradSolVAR[ivar][jvar];
	      //resRhs  += SolVAR[jvar]*gradSolVAR[ivar][jvar]-0*IRe*NablaSolVAR[ivar][jvar];
	      supgPhi += (SolVAR[jvar]*gradphi2[i*dim+jvar])* supgTau; 
	    }
	    //resRhs += GradSolP[ivar];
	    F[SolPdeIndex[ivar]][i]+= ( -IRe*Lap_rhs-NavierStokes*Adv_rhs[ivar]*(phi2[i]+supgPhi)
					+SolVAR[dim]*gradphi2[i*dim+ivar])*Weight2;
	  }
	  //END RESIDUALS A block ===========================
	  
	  if(assembe_matrix){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve2; j++) {
	      double Lap=0;
	      
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		// Laplacian
		Lap  += (gradphi2[i*dim+ivar]*gradphi2[j*dim+ivar])*Weight2;
	      }
	      for(unsigned ivar=0; ivar<dim; ivar++) {    
		B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][i*nve2+j] += IRe*Lap;
		  		  		  
		  for(unsigned jvar=0; jvar<dim; jvar++) {
		    B[SolPdeIndex[ivar]][SolPdeIndex[ivar]][i*nve2+j] += (phi2[i]+supgPhi)*SolVAR[jvar]*gradphi2[j*dim+jvar]*Weight2;
		    B[SolPdeIndex[ivar]][SolPdeIndex[jvar]][i*nve2+j] += (phi2[i]+supgPhi)*phi2[j]*gradSolVAR[ivar][jvar]*Weight2;
		    B[SolPdeIndex[ivar]][SolPdeIndex[jvar]][i*nve2+j] +=  phi2[j]*gradphi2[i*dim+jvar]*supgTau*Adv_rhs[ivar]*Weight2;
		  }
		  
		
	      }
  	    } //end phij loop
	    
	    // *** phi1_j loop ***
	    for(unsigned j=0; j<nve1; j++){
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndex[ivar]][SolPdeIndex[dim]][i*nve1+j] -= gradphi2[i*dim+ivar]*phi1[j]*Weight2;
	      }
	    } //end phi1_j loop
	  } // endif assembe_matrix
	} //end phii loop
  

	// *** phi1_i loop ***
	for(unsigned i=0; i<nve1; i++){
	  //BEGIN RESIDUALS B block ===========================
	  double div = 0;
	  for(unsigned ivar=0; ivar<dim; ivar++) {
	    div += gradSolVAR[ivar][ivar];
	  }
	  F[SolPdeIndex[dim]][i]+= (phi1[i]*div)*Weight2;// + 
	                             //+ 0*hk*hk*(1./0.001)*alpha*(GradSolP[0]*gradphi2[i*dim + 0] + GradSolP[1]*gradphi2[i*dim + 1]) )*Weight2;
				     
           
	  //END RESIDUALS  B block ===========================
	  
	  if(assembe_matrix){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve2; j++) {
	      for(unsigned ivar=0; ivar<dim; ivar++) {
		B[SolPdeIndex[dim]][SolPdeIndex[ivar]][i*nve2+j]-= phi1[i]*gradphi2[j*dim+ivar]*Weight2;
	      }
	    }  //end phij loop
	  } // endif assembe_matrix
	}  //end phi1_i loop
	
	if(assembe_matrix * penalty){  //block nve1 nve1
	  // *** phi_i loop ***
	  for(unsigned i=0; i<nve1; i++){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve1; j++){
	      //B[SolPdeIndex[dim]][SolPdeIndex[dim]][i*nve1+j]-= ILambda*phi1[i]*phi1[j]*Weight2;
	      for(unsigned ivar=0; ivar<dim; ivar++) {
	        B[SolPdeIndex[dim]][SolPdeIndex[dim]][i*nve1+j] -= 0.;//*hk*hk*(1./0.001)*alpha*(gradphi2[i*dim + ivar]*gradphi2[j*dim + ivar])*Weight2; 
	      }
	    }
	  }
	}   //end if penalty
      }  // end gauss point loop
      
      //--------------------------------------------------------------------------------------------------------
      // Boundary Integral --> to be added
      //number of faces for each type of element
//       if (igrid==gridn || !myel->GetRefinedElementIndex(kel) ) {
//      
// 	unsigned nfaces = myel->GetElementFaceNumber(kel);
// 
// 	// loop on faces
// 	for(unsigned jface=0;jface<nfaces;jface++){ 
// 	  
// 	  // look for boundary faces
// 	  if(myel->GetFaceElementIndex(kel,jface)<0){
// 	    for(unsigned ivar=0; ivar<dim; ivar++) {
// 	      ml_prob.ComputeBdIntegral(pdename, &Solname[ivar][0], kel, jface, level, ivar);
// 	    }
// 	  }
// 	}	
//       }
      //--------------------------------------------------------------------------------------------------------
    } // endif single element not refined or fine grid loop
    //--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the Global Matrix/Vector
    for(unsigned ivar=0; ivar<dim; ivar++) {
      myRES->add_vector_blocked(F[SolPdeIndex[ivar]],KK_dof[ivar]);
      if(assembe_matrix){
	myKK->add_matrix_blocked(B[SolPdeIndex[ivar]][SolPdeIndex[ivar]],KK_dof[ivar],KK_dof[ivar]);  
	myKK->add_matrix_blocked(B[SolPdeIndex[ivar]][SolPdeIndex[dim]],KK_dof[ivar],KK_dof[dim]);
	myKK->add_matrix_blocked(B[SolPdeIndex[dim]][SolPdeIndex[ivar]],KK_dof[dim],KK_dof[ivar]);
	if(nwtn_alg==2){
	  for(unsigned jvar=1; jvar<dim; jvar++) {
	    myKK->add_matrix_blocked(B[SolPdeIndex[ivar]][SolPdeIndex[(ivar+jvar)%dim]],KK_dof[ivar],KK_dof[(ivar+jvar)%dim]);  
	  }
	}
      }
    }
    //Penalty
    if(assembe_matrix*penalty) myKK->add_matrix_blocked(B[SolPdeIndex[dim]][SolPdeIndex[dim]],KK_dof[dim],KK_dof[dim]);
    myRES->add_vector_blocked(F[SolPdeIndex[dim]],KK_dof[dim]);
    //--------------------------------------------------------------------------------------------------------  
  } //end list of elements loop for each subdomain
  
  
  if(assembe_matrix) myKK->close();
  myRES->close();
  // ***************** END ASSEMBLY *******************
}*/

//------------------------------------------------------------------------------------------------------------
void AssembleMatrixResT(MultiLevelProblem &ml_prob, unsigned level, const unsigned &gridn, const bool &assembe_matrix){
  
  //pointers and references
  Solution*      mysolution	       = ml_prob._ml_sol->GetSolutionLevel(level);
  LinearImplicitSystem& mylin_impl_sys = ml_prob.get_system<LinearImplicitSystem>("Temperature");
  LinearEquationSolver*  mylsyspde     = mylin_impl_sys._LinSolver[level];   
  mesh*          mymsh		       = ml_prob._ml_msh->GetLevel(level);
  elem*          myel		       = mymsh->el;
  SparseMatrix*  myKK		       = mylsyspde->_KK;
  NumericVector* myRES		       = mylsyspde->_RES;
  MultiLevelSolution* ml_sol           = ml_prob._ml_sol;
  
  
  //data
  const unsigned	dim	= mymsh->GetDimension();
  unsigned 		nel	= mymsh->GetElementNumber();
  unsigned 		igrid	= mymsh->GetGridNumber();
  unsigned 		iproc	= mymsh->processor_id();
  double		IPe	= 1./(ml_prob.parameters.get<Fluid>("Fluid").get_Peclet_number());  
  
  //solution variable
  unsigned SolIndex;  
  unsigned SolPdeIndex;
  SolIndex=ml_sol->GetIndex("T");
  SolPdeIndex=mylin_impl_sys.GetSolPdeIndex("T");
  //solution order
  unsigned order_ind = ml_sol->GetSolutionType(SolIndex);
  unsigned end_ind   = mymsh->GetEndIndex(order_ind);
  
  //coordinates
  vector< vector < double> > coordinates(dim); 
  //const char coordinate_name[3][2] = {"X","Y","Z"};
  //vector < unsigned > coordinate_Index(dim);
//   for(unsigned ivar=0; ivar<dim; ivar++) {
//     coordinate_Index[ivar]=ivar;//ml_prob.GetIndex(coordinate_name[ivar]);
//   }
  
  // declare 
  vector< int > metis_node;
  vector< int > KK_dof;
  vector <double> phi;
  vector <double> gradphi;
  vector <double> nablaphi;
    
  double weight;
  vector< double > F;
  vector< double > B;
 
  // reserve 
  const unsigned max_size = static_cast< unsigned > (ceil(pow(3,dim)));
  metis_node.reserve(max_size);
  KK_dof.reserve(max_size);
  for(int i=0;i<dim;i++) 
    coordinates[i].reserve(max_size);
  phi.reserve(max_size);
  gradphi.reserve(max_size*dim);
  nablaphi.reserve(max_size*(3*(dim-1)));
  F.reserve(max_size);
  B.reserve(max_size*max_size);
  
  // Set to zeto all the entries of the Global Matrix
  if(assembe_matrix) myKK->zero();
  
  // *** element loop ***
 
  for (int iel=mymsh->IS_Mts2Gmt_elem_offset[iproc]; iel < mymsh->IS_Mts2Gmt_elem_offset[iproc+1]; iel++) {

    unsigned kel = mymsh->IS_Mts2Gmt_elem[iel];
    short unsigned kelt=myel->GetElementType(kel);
    unsigned nve=myel->GetElementDofNumber(kel,end_ind);
    
    // resize
    metis_node.resize(nve);
    KK_dof.resize(nve);
    phi.resize(nve);
    gradphi.resize(nve*dim);
    nablaphi.resize(nve*(3*(dim-1)));	
    for(int i=0;i<dim;i++){
      coordinates[i].resize(nve);
    }
    
    // set to zero all the entries of the FE matrices
    F.resize(nve);
    memset(&F[0],0,nve*sizeof(double));
    if(assembe_matrix){
      B.resize(nve*nve);
      memset(&B[0],0,nve*nve*sizeof(double));
    }
    
    // get local to global mappings
    for( unsigned i=0;i<nve;i++){
      unsigned inode=myel->GetElementVertexIndex(kel,i)-1u;	
      unsigned inode_coord_metis=mymsh->GetMetisDof(inode,2);
      metis_node[i]=mymsh->GetMetisDof(inode,order_ind);
      for(unsigned ivar=0; ivar<dim; ivar++) {
	coordinates[ivar][i]=(*mymsh->_coordinate->_Sol[ivar])(inode_coord_metis);
      }
      KK_dof[i]=mylsyspde->GetKKDof(SolIndex,SolPdeIndex,inode);
    }
        
    if(igrid==gridn || !myel->GetRefinedElementIndex(kel)) {
      // *** Gauss poit loop ***
      for(unsigned ig=0;ig < ml_prob._ml_msh->_finiteElement[kelt][order_ind]->GetGaussPointNumber(); ig++) {
	// *** get Jacobian and test function and test function derivatives ***
	ml_prob._ml_msh->_finiteElement[kelt][order_ind]->Jacobian(coordinates,ig,weight,phi,gradphi,nablaphi);
	//Temperature and velocity current solution
	double SolT=0;
	vector < double > gradSolT(dim,0.);
	for(unsigned ivar=0; ivar<dim; ivar++){
	  gradSolT[ivar]=0; 
	}
	vector < double > SolU(dim,0.);
	vector < unsigned > SolIndexU(dim);
	SolIndexU[0]=ml_sol->GetIndex("U");
	SolIndexU[1]=ml_sol->GetIndex("V");
	if(dim==3) SolIndexU[2]=ml_sol->GetIndex("W");
	  	  
	unsigned SolType=ml_sol->GetSolutionType("T");
	for(unsigned i=0; i<nve; i++) {
	  double soli = (*mysolution->_Sol[SolIndex])(metis_node[i]);
	  SolT+=phi[i]*soli;
	  for(unsigned ivar=0; ivar<dim; ivar++) gradSolT[ivar] += gradphi[i*dim+ivar]*soli; 
	  for(int j=0;j<dim;j++)  {
	    SolU[j]+=phi[i]*(*mysolution->_Sol[SolIndexU[j]])(metis_node[i]);
	  }
	}
	// *** phi_i loop ***
	for(unsigned i=0; i<nve; i++){
	  //BEGIN RESIDUALS A block ===========================
	  double Adv_rhs=0;
	  double Lap_rhs=0;
	  for(unsigned ivar=0; ivar<dim; ivar++) {
	    Lap_rhs += gradphi[i*dim+ivar]*gradSolT[ivar];
	    Adv_rhs += SolU[ivar]*gradSolT[ivar];
	  }
	  F[i]+= (-IPe*Lap_rhs-Adv_rhs*phi[i])*weight; 		    
	  //END RESIDUALS A block ===========================
	  if(assembe_matrix){
	    // *** phi_j loop ***
	    for(unsigned j=0; j<nve; j++) {
	      double Lap=0;
	      double Adv1=0;
	      for(unsigned ivar=0; ivar<dim; ivar++) {
	      // Laplacian
		Lap  += gradphi[i*dim+ivar]*gradphi[j*dim+ivar]*weight;
		// advection term I
		Adv1 += SolU[ivar]*gradphi[j*dim+ivar]*phi[i]*weight;
	      }
	      B[i*nve+j] += IPe*Lap + Adv1;
	    } // end phij loop
	  } // end phii loop
	} // endif assembe_matrix
      } // end gauss point loop
    } // endif single element not refined or fine grid loop
    //--------------------------------------------------------------------------------------------------------
    //Sum the local matrices/vectors into the global Matrix/vector
      
    myRES->add_vector_blocked(F,KK_dof);
    if(assembe_matrix) myKK->add_matrix_blocked(B,KK_dof,KK_dof);  
  } //end list of elements loop for each subdomain
    
  myRES->close();
  if(assembe_matrix) myKK->close();
  
   // ***************** END ASSEMBLY *******************
  
}

