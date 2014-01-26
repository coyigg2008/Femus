//#include "SolverlibConf.hpp"
#include "FEMTTUConfig.h"

#ifdef HAVE_PETSC

// Local Includes
#include "PetscMacro.hpp"
#include "PetscLinearSolver.hpp"
#include "PetscPreconditioner.hpp"
#include "PetscVector.hpp"
#include "PetscMatrix.hpp"
#include "PetscRectangularMatrix.hpp"

using namespace std;

// ==========================================================
extern "C" {
#if PETSC_VERSION_LESS_THAN(2,2,1)
  typedef int PetscErrorCode;
  typedef int PetscInt;
#endif

#if PETSC_VERSION_LESS_THAN(3,0,1) && PETSC_VERSION_RELEASE
  // ------------------------------------------------------------------
  PetscErrorCode __libmesh_petsc_preconditioner_setup(void * ctx) {
    Preconditioner * preconditioner = static_cast<Preconditioner*>(ctx);
    preconditioner->init();
    return 0;
  }
// ------------------------------------------------------------------------------
  PetscErrorCode __libmesh_petsc_preconditioner_apply(void *ctx, Vec x, Vec y)  {
    Preconditioner * preconditioner = static_cast<Preconditioner*>(ctx);
    PetscVector x_vec(x);
    PetscVector y_vec(y);
    preconditioner->apply(x_vec,y_vec);
    return 0;
  }
#else
// ----------------------------------------------------------------
  PetscErrorCode __libmesh_petsc_preconditioner_setup(PC pc) {
    void *ctx;
    PetscErrorCode ierr = PCShellGetContext(pc,&ctx);
    CHKERRQ(ierr);
    Preconditioner * preconditioner = static_cast<Preconditioner*>(ctx);
    preconditioner->init();
    return 0;
  }
// --------------------------------------------------------------------------
  PetscErrorCode __libmesh_petsc_preconditioner_apply(PC pc, Vec x, Vec y) {
    void *ctx;
    PetscErrorCode ierr = PCShellGetContext(pc,&ctx);
    CHKERRQ(ierr);
    Preconditioner * preconditioner = static_cast<Preconditioner*>(ctx);
    PetscVector x_vec(x);
    PetscVector y_vec(y);
    preconditioner->apply(x_vec,y_vec);
    return 0;
  }
#endif
} // end extern "C
// ================================================


// ==============================================
// ----------------------- functions ------
// ==============================================

void PetscLinearSolver::set_tolerances(const double rtol, const double atol,
                                        const double divtol, const unsigned maxits) {

  _rtol   = static_cast<PetscReal>(rtol);
  _abstol = static_cast<PetscReal>(atol);
  _dtol   = static_cast<PetscReal>(divtol);
  _maxits = static_cast<PetscInt>(maxits);

};


void PetscLinearSolver::set_schur_tolerances(const double rtol, const double atol,
    const double divtol, const unsigned maxits) {

  _rtol_   = static_cast<PetscReal>(rtol);
  _abstol_ = static_cast<PetscReal>(atol);
  _dtol_   = static_cast<PetscReal>(divtol);
  _maxits_ = static_cast<PetscInt>(maxits);

};


int PetscLinearSolver::Vanka_Smoother(const vector <unsigned> &MGIndex,const vector <unsigned> &VankaIndex,
                                       const short unsigned &NSchurVar,const bool &Schur) {

  SearchTime=0, AssemblyTime=0, SolveTime0=0, SolveTime1=0, SolveTime2=0, UpdateTime=0;

  its_A=0; its_C=0; its=0;

  unsigned nvt=KKIndex[KKIndex.size()-1u];
  unsigned nel=_msh->GetElementNumber();
  int grid=_msh->GetGridNumber();
  bool FastVankaBlock=(_SolType[MGIndex[VankaIndex[VankaIndex.size()-NSchurVar]]]<3)?false:true;
  
  if (indexa.size()<nvt || indexc.size()<nel) {
    indexa.resize(nvt);
    indexb.resize(nvt);
    indexc.resize(nel);
    indexd.resize(nel);
  }
  
  long unsigned IndexaSize=nvt;
  long unsigned IndexbSize=nvt;
  long unsigned IndexcSize=nel;
  long unsigned IndexdSize=nel;

  if (indexai.size()<IndexaSize) {
    indexai.resize(IndexaSize);
    indexbi.resize(IndexbSize);
    indexci.resize(IndexcSize);
    indexdi.resize(IndexdSize);
  }

  for (unsigned i=0; i<nvt; i++) {
    indexa[i]=IndexaSize;
    indexb[i]=IndexbSize;
  }
  for (unsigned i=0; i<nel; i++) {
    indexc[i]=IndexcSize;
    indexd[i]=IndexdSize;
  }

  /// *** Start Vanka Block ***
  for(int isdom=0; isdom<_msh->nsubdom; isdom++) {
//    
//    unsigned isdom = mapsubdomain[iasdom];
   
//  for (int gel=0; gel < static_cast <int> (nel); gel += _num_elem_vanka_block) {
  for (int gel=_msh->IS_Mts2Gmt_elem_offset[isdom]; gel < _msh->IS_Mts2Gmt_elem_offset[isdom+1]; gel += _num_elem_vanka_block) {
    // ***************** NODE/ELEMENT SERCH *******************
    start_time=clock();
    Asize=0;
    counterb=0;
    Csize=0;
    Dsize=0;
    PDsize=0;
    
//     cout << gel << endl;
//     cout << "isdom " << isdom << "  " << IS_Mts2Gmt_elem_offset[isdom] << "  " << IS_Mts2Gmt_elem_offset[isdom+1] << endl;
//      cout << "+++++++++++" << endl;
//     cout << "metis offset: " << MetisOffset[2][isdom+1] << MetisOffset[4][isdom+1] << endl;
    
//      for (int iel=gel; iel<gel+_num_elem_vanka_block && iel< nel; iel++) {
     for (int iel_mts=gel; iel_mts<gel+_num_elem_vanka_block && iel_mts< _msh->IS_Mts2Gmt_elem_offset[isdom+1]; iel_mts++) {
      
//        cout << "gel: " << gel << "  iel: " << iel_mts << endl; 
      // IS_Mts2Gmt_elem from Metis dof to Gmt element dof     
//       cout << "iel metis: " << iel << "  iel gambit " << IS_Mts2Gmt_elem[iel] << endl;
      
      unsigned iel = _msh->IS_Mts2Gmt_elem[iel_mts];     
//      cout << "iel metis: " << iel_mts << "  iel gambit " << iel << endl;
      
      for (unsigned i=0; i<_msh->el->GetElementDofNumber(iel,0); i++) {
        unsigned inode=_msh->el->GetElementVertexIndex(iel,i)-1u;
        unsigned nvei=_msh->el->GetVertexElementNumber(inode);
        const unsigned *pt_jel=_msh->el->GetVertexElementAddress(inode,0);
        for (unsigned j=0; j<nvei*(!FastVankaBlock)+FastVankaBlock; j++) {
          unsigned jel=(!FastVankaBlock)?*(pt_jel++)-1u:iel;
          //add elements for velocity to be solved
          if (indexc[jel]==IndexcSize) {
            indexci[Csize]=jel;
            indexc[jel]=Csize++;
	    //----------------------------------------------------------------------------------
            ///add velocity nodes to be solved
            for (unsigned iind=0; iind<VankaIndex.size()-NSchurVar; iind++) {
              unsigned indexSol=VankaIndex[iind];
              const unsigned *pt_un=_msh->el->GetElementVertexAddress(jel,0);
              unsigned nvej=_msh->el->GetElementDofNumber(jel,_msh->_END_IND[_SolType[MGIndex[indexSol]]]);
              for (unsigned jj=0; jj<nvej; jj++) {
                unsigned jnode=(_SolType[MGIndex[indexSol]]<3)?(*(pt_un++)-1u):(jel+jj*nel);
		unsigned kkdof=GetKKDof(MGIndex[indexSol], indexSol, jnode);
		unsigned jnode_Metis = _msh->GetMetisDof(jnode,_SolType[MGIndex[indexSol]]);
		if (indexa[kkdof]==IndexaSize && 1.1 <(*(*_Bdc)[MGIndex[indexSol]])(jnode_Metis) ) {
                  indexai[Asize]=kkdof;
                  indexa[kkdof]=Asize++;
		}
              }
            }
            //-----------------------------------------------------------------------------------
            for (unsigned jj=0; jj<_msh->el->GetElementDofNumber(jel,0); jj++) {
              unsigned jnode=_msh->el->GetElementVertexIndex(jel,jj)-1u;
              unsigned nvej=_msh->el->GetVertexElementNumber(jnode);
              const unsigned *pt_kel=_msh->el->GetVertexElementAddress(jnode,0);
              for (unsigned k=0; k<nvej; k++) {
                unsigned kel=*(pt_kel++)-1u;
                ///add all variables to be updated
                if (indexd[kel]==IndexdSize) {
                  indexdi[Dsize]=kel;
                  indexd[kel]=Dsize++;
                 for (unsigned int indexSol=0; indexSol<KKIndex.size()-1u; indexSol++) {
                    const unsigned *pt_un=_msh->el->GetElementVertexAddress(kel,0);
                    unsigned nvek=_msh->el->GetElementDofNumber(kel,_msh->_END_IND[_SolType[MGIndex[indexSol]]]);
                    for (unsigned kk=0; kk<nvek; kk++) {
                      unsigned knode=(_SolType[MGIndex[indexSol]]<3)?(*(pt_un++)-1u):(kel+kk*nel);
		      unsigned kkdof=GetKKDof(MGIndex[indexSol], indexSol, knode);
		      unsigned knode_Metis = _msh->GetMetisDof(knode,_SolType[MGIndex[indexSol]]);
		       if (indexb[kkdof]==IndexbSize && 0.1<(*(*_Bdc)[MGIndex[indexSol]])(knode_Metis)) {
                        indexbi[counterb]=kkdof;
                        indexb[kkdof]=counterb++;
                      }
                    }
                  }
                }
              }
            } //------------------------------------------------------------------------
          }
        }
      }
      //-----------------------------------------------------------------------------------------
      ///Add pressure nodes to be solved
      {
        for (unsigned iind=VankaIndex.size()-NSchurVar; iind<VankaIndex.size(); iind++) {
          unsigned indexSol=VankaIndex[iind];
          const unsigned *pt_un=_msh->el->GetElementVertexAddress(iel,0);
          unsigned nvei=_msh->el->GetElementDofNumber(iel,_msh->_END_IND[_SolType[MGIndex[indexSol]]]);
          for (unsigned ii=0; ii<nvei; ii++) {
            unsigned inode=(_SolType[MGIndex[indexSol]]<3)?(*(pt_un++)-1u):(iel+ii*nel);
	    unsigned kkdof=GetKKDof(MGIndex[indexSol], indexSol, inode);
	    unsigned inode_Metis = _msh->GetMetisDof(inode,_SolType[MGIndex[indexSol]]);
	    if (indexa[kkdof]==IndexaSize && 1.1<(*(*_Bdc)[MGIndex[indexSol]])(inode_Metis) ) {
              indexai[Asize]=kkdof;
              indexa[kkdof]=Asize++;
	      PDsize++;
	    }
          }
        }
      }
      //-----------------------------------------------------------------------------------------
    }

    PBsize=Asize;
    PCsize=PBsize-PDsize;
    for (PetscInt i=0; i<counterb; i++) {
      unsigned jnode=indexbi[i];
      if (indexa[jnode]==IndexaSize) {
        indexai[Asize]=jnode;
        indexa[jnode]=Asize++;
      }
      // *** reinitialize indexb
      indexb[jnode]=IndexbSize;
    }
    PAmBsize=Asize-PBsize;
    
    // *** re-initialize indeces(a,c,d)
    for (PetscInt i=0; i<Asize; i++) {
      indexa[indexai[i]]=IndexaSize;
    }
    for (PetscInt i=0; i<Csize; i++) {
      indexc[indexci[i]]=IndexcSize;
    }
    for (PetscInt i=0; i<Dsize; i++) {
      indexd[indexdi[i]]=IndexdSize;
    }

    end_time=clock();
    SearchTime+=(end_time-start_time);
    start_time=clock();
    /// ***************** END NODE/ELEMENT SEARCH *******************

    /// ***************** ASSEMBLY *******************

    // generate IS
    ierr=ISCreateGeneral(PETSC_COMM_SELF,PBsize,&indexai[0],PETSC_USE_POINTER,&isPA);
    CHKERRQ(ierr);
    ierr=ISSort(isPA); 
    CHKERRQ(ierr);
    
    // Solution Vector Pw
    ierr = VecCreateSeq(PETSC_COMM_SELF,PBsize,&Pw);
    CHKERRQ(ierr);
    ierr = VecSetFromOptions(Pw);
    CHKERRQ(ierr);

    // RHS Pr
    ierr = VecDuplicate(Pw,&Pr);
    CHKERRQ(ierr);

    //Matrix PA (PA is the Total matrix over a patch)
    ierr = MatGetSubMatrix(KK,isPA,isPA,MAT_INITIAL_MATRIX,&PA);
    CHKERRQ(ierr);

    if (!Schur || grid==0) {
        
      ISCreateStride(PETSC_COMM_SELF,PBsize,0,1,&isPB);
      PetscScalar *y=new PetscScalar [PBsize];

      ierr = ISGetIndices(isPA,&ind);
      CHKERRQ(ierr);
      ierr = ISGetIndices(isPB,&ind2);
      CHKERRQ(ierr);
      ierr = VecGetValues(RES,PBsize,ind,y);
      ierr = VecSetValues(Pr,PBsize,ind2,y,INSERT_VALUES);
      CHKERRQ(ierr);
      ierr = VecAssemblyBegin(Pr);
      CHKERRQ(ierr);
      ierr = VecAssemblyEnd(Pr);
      CHKERRQ(ierr);
      ierr = ISRestoreIndices(isPA,&ind);
      CHKERRQ(ierr);
      ierr = ISRestoreIndices(isPB,&ind2);
      CHKERRQ(ierr);
      ierr = ISDestroy(&isPB);
      CHKERRQ(ierr);  
      
      delete [] y;

      end_time=clock();
      AssemblyTime+=(end_time-start_time);
      /// ***************** END ASSEMBLY ******************
      

      /// ***************** START NON SCHUR COMPLEMENT SOLVER ******************
      start_time=clock(); // START SOLVE 0 TIME
//        PetscViewer    viewer;
//        ierr=PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,600,600,&viewer);CHKERRQ(ierr);
// //        ierr= MatView(PA,viewer);CHKERRQ(ierr);
//        ierr= VecView(Pr,viewer);CHKERRQ(ierr);
//        double ff;
//        std::cin>>ff;
//        exit(0);

      // Initialize the data structure with the matrix PA
      this->init(PA,false,false);
      
//       ierr=VecZeroEntries(Pw);
//       CHKERRQ(ierr);
      
      // Solve the linear system
      ierr = KSPSolve(_ksp,Pr,Pw);   //CHKERRQ(ierr);

      // print information on ksp solver
      //ierr = KSPView(_ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

      ierr = KSPGetIterationNumber(_ksp,&its);
      its_A += its;

      end_time=clock();
      SolveTime0+=(end_time-start_time);

      /// ***************** END NON SCHUR COMPLEMENT SOLVER ******************
    }

    else {

      /// ***************** START SCHUR COMPLEMENT SOLVER ******************
   
      // residual w (Velocity only)
      ierr = VecCreateSeq(PETSC_COMM_SELF,PCsize,&w);
      CHKERRQ(ierr);
      ierr = VecSetFromOptions(w);
      CHKERRQ(ierr);
      // RHS f (Momentum only)
      ierr = VecDuplicate(w,&f);
      CHKERRQ(ierr);

      // residual z (Pressure only)
      ierr = VecCreateSeq(PETSC_COMM_SELF,PDsize,&z);
      CHKERRQ(ierr);
      ierr = VecSetFromOptions(z);
      CHKERRQ(ierr);
      // RHS g (Continuity only);
      ierr = VecDuplicate(z,&g);
      CHKERRQ(ierr);
      
      ISCreateStride(PETSC_COMM_SELF,PBsize,0,1,&isPB);
      ierr = ISGetIndices(isPB,&ind2);
      
      //da sistemare (TODO I think that allocation is not necessary)
      PetscScalar *y=new PetscScalar [PBsize];
      
      ierr = VecGetValues(RES,PBsize,&indexai[0],y);
      ierr = VecSetValues(f,PCsize,ind2,y,INSERT_VALUES);
      CHKERRQ(ierr);
      ierr = VecAssemblyBegin(f);
      CHKERRQ(ierr);
      ierr = VecAssemblyEnd(f);
      CHKERRQ(ierr);
      ierr = VecSetValues(g,PDsize,ind2,&y[PCsize],INSERT_VALUES);
      CHKERRQ(ierr);
      ierr = VecAssemblyBegin(g);
      CHKERRQ(ierr);
      ierr = VecAssemblyEnd(g);
      CHKERRQ(ierr);

      delete [] y;

      ierr=ISCreateGeneral(PETSC_COMM_SELF,PCsize,ind2,PETSC_USE_POINTER,&isCl);
      CHKERRQ(ierr);
      ierr=ISCreateGeneral(PETSC_COMM_SELF,PDsize,&ind2[PCsize],PETSC_USE_POINTER,&isDl);
      CHKERRQ(ierr);

      /* linear system matrix Schur Complement*/
      ierr = MatGetSubMatrix(PA,isCl,isCl,MAT_INITIAL_MATRIX,&A);
      CHKERRQ(ierr);
      ierr = MatGetSubMatrix(PA,isDl,isCl,MAT_INITIAL_MATRIX,&B);
      CHKERRQ(ierr);
      
      if (GetMatrixProperties()) {
        ierr = MatCreateTranspose(B,&Bt);
        CHKERRQ(ierr);
      } else {
        // In the Monolithic Fluid-Structure formulation Bt is not the transpose of B
        ierr = MatGetSubMatrix(PA,isCl,isDl,MAT_INITIAL_MATRIX,&Bt);
        CHKERRQ(ierr);
      }

      //C=D-B A^{-1} B^t
      if (!GetStabilization()) {
        ierr = MatCreateSchurComplement(A,A,Bt,B,PETSC_NULL,&C);
        CHKERRQ(ierr);
      } else {
        ierr = MatGetSubMatrix(PA,isDl,isDl,MAT_INITIAL_MATRIX,&D);
        CHKERRQ(ierr);
        ierr = MatCreateSchurComplement(A,A,Bt,B,D,&C);
        CHKERRQ(ierr);
      }

      ierr = ISDestroy(&isCl);
      CHKERRQ(ierr);
      ierr = ISDestroy(&isDl);
      CHKERRQ(ierr);
      
      end_time=clock();
      AssemblyTime+=(end_time-start_time);

      /// ***************** END ASSEMBLY ******************

      /// ***************** START SCHUR COMPLEMENT SOLVER ******************

      start_time=clock(); // START SOLVE 0 TIME

      //init Schur Complement
      this->init_schur(C);

      //init _Ksp (GMRESM - no precond)
      this->init(C,true,Schur);
       
/*      PetscViewer    viewer;
      ierr=PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,600,600,&viewer);CHKERRQ(ierr);
      ierr= MatView(A,viewer);CHKERRQ(ierr);
      double ff;
      std::cin>>ff;  */    
      

      /// Prediction step
//       ierr=VecZeroEntries(w);
//       CHKERRQ(ierr);
//       
      ierr = KSPSolve(_ksp2,f,w);
      CHKERRQ(ierr);   // solve A^{-1} f -> w (temporary)

      ierr = KSPGetIterationNumber(_ksp2,&its);
      its_A += its;

      end_time=clock();  // END SOLVE 0 TIME
      SolveTime0+=(end_time-start_time);
      start_time=clock(); // START SOLVE 1 TIME

      /// Projection Step
      // Computation of the Rhs for the Projection Step
      ierr = VecScale (w, -1.);
      CHKERRQ(ierr);
      ierr = MatMultAdd(B,w,g,g);
      CHKERRQ(ierr); // evaluate g - B A^{-1} f -> g (temporary)

//       ierr=VecZeroEntries(z);
//       CHKERRQ(ierr);
      
      ierr = KSPSolve(_ksp,g,z);
      CHKERRQ(ierr); // solve z = C^{-1} (g - B A^{-1} f)

      ierr = KSPGetIterationNumber(_ksp,&its);
      its_C += its;

      // print information on ksp solver
      // ierr = KSPView(_ksp2,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

      end_time=clock(); // END SOLVE 1 TIME
      SolveTime1+=(end_time-start_time);
      start_time=clock(); // START SOLVE 2 TIME

      /// Correction Step
      // Computation of the Rhs for the Correction Step
      ierr =VecCopy(z,g);
      CHKERRQ(ierr);
      ierr = VecScale(g, -1.);
      CHKERRQ(ierr);
      ierr = MatMultAdd(Bt,g,f,f);
      CHKERRQ(ierr);  // evaluate f - Bt z -> f

      ierr = KSPSolve(_ksp2,f,w);
      CHKERRQ(ierr); // solve w=A^{-1}(f-Bt z)

      PetscScalar *Z[1];
      ierr = VecGetArray(w,Z);
      CHKERRQ(ierr);
      ierr = VecSetValues(Pw,PCsize,ind2,Z[0],INSERT_VALUES);
      CHKERRQ(ierr);
      ierr = VecRestoreArray(w,Z);
      CHKERRQ(ierr);

      ierr = VecGetArray(z,Z);
      CHKERRQ(ierr);
      ierr = VecSetValues(Pw,PDsize,&ind2[PCsize],Z[0],INSERT_VALUES);
      CHKERRQ(ierr);
      ierr = VecRestoreArray(z,Z);
      CHKERRQ(ierr);

      ierr = VecAssemblyBegin(Pw);
      CHKERRQ(ierr);
      ierr = VecAssemblyEnd(Pw);
      CHKERRQ(ierr);
      
      ierr=ISRestoreIndices(isPB,&ind2); 
      CHKERRQ(ierr);

      ierr = ISDestroy(&isPB);
      CHKERRQ(ierr);
      end_time=clock(); // END SOLVE 2 TIME
      SolveTime2+=(end_time-start_time);


    }
    /// ***************** END VANKA SCHUR COMPLEMENT SOLVER ******************

    /// ***************** START UPDATE ******************

    start_time=clock();

    //update Residual for PA
    ierr = MatMult(PA,Pw,Pr);
    CHKERRQ(ierr);

    ierr = VecScale(Pr, -1.);
    CHKERRQ(ierr);
    PetscScalar *R[1];
    ierr = VecGetArray(Pr,R);
    CHKERRQ(ierr);

    PetscScalar *W[1];
    ierr = VecGetArray(Pw,W);
    CHKERRQ(ierr);

    ierr = ISGetIndices(isPA,&ind);
    CHKERRQ(ierr);

    ierr=VecSetValues(RES,PBsize,ind,R[0],ADD_VALUES);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(RES);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(RES);
    CHKERRQ(ierr);

    ierr=VecSetValues(EPS,PBsize,ind,W[0],ADD_VALUES);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(EPS);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(EPS);
    CHKERRQ(ierr);

    ierr = ISRestoreIndices(isPA,&ind);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(Pr,R);
    CHKERRQ(ierr);

    if (PAmBsize) {
      ierr=ISCreateGeneral(PETSC_COMM_SELF,PAmBsize,&indexai[PBsize],PETSC_USE_POINTER,&isB);
      CHKERRQ(ierr);
      ierr=ISSort(isB);
      CHKERRQ(ierr);

      ierr = VecCreateSeq(PETSC_COMM_SELF,PAmBsize,&Ps);
      CHKERRQ(ierr);

      ierr = VecSetFromOptions(Ps);
      CHKERRQ(ierr);

      ierr = MatGetSubMatrix(KK,isB,isPA,MAT_INITIAL_MATRIX,&PB);
      CHKERRQ(ierr);
      ierr = MatMult(PB,Pw,Ps);
      CHKERRQ(ierr);

      ierr = VecScale (Ps, -1.);
      CHKERRQ(ierr);
      PetscScalar *S[1];
      ierr = VecGetArray(Ps,S);
      CHKERRQ(ierr);
      ierr = ISGetIndices(isB,&ind);
      CHKERRQ(ierr);

      ierr=VecSetValues(RES,PAmBsize,ind,S[0],ADD_VALUES);
      CHKERRQ(ierr);
      ierr = VecAssemblyBegin(RES);
      CHKERRQ(ierr);
      ierr = VecAssemblyEnd(RES);
      CHKERRQ(ierr);

      ierr = ISRestoreIndices(isB,&ind);
      CHKERRQ(ierr);
      ierr = VecRestoreArray(Ps,S);
      CHKERRQ(ierr);
      ierr = VecDestroy(&Ps);
      CHKERRQ(ierr);
      ierr = MatDestroy(&PB);
      CHKERRQ(ierr);
      ierr = ISDestroy(&isB);
      CHKERRQ(ierr);
    }

    /// ***************** START CLEANING ******************
    
    ierr = VecRestoreArray(Pw,W);
    CHKERRQ(ierr);

    ierr = VecDestroy(&Pw);
    CHKERRQ(ierr);
    ierr = VecDestroy(&Pr);
    CHKERRQ(ierr);
    ierr = MatDestroy(&PA);
    CHKERRQ(ierr);
    ierr = ISDestroy(&isPA);
    CHKERRQ(ierr);

    ierr = KSPDestroy(&_ksp);
    CHKERRQ(ierr);
    if (Schur && grid > 0) {
      ierr = MatDestroy(&A);
      CHKERRQ(ierr);
      ierr = MatDestroy(&B);
      CHKERRQ(ierr);
      ierr = MatDestroy(&Bt);
      CHKERRQ(ierr);
      ierr = MatDestroy(&C);
      CHKERRQ(ierr);
      if (GetStabilization()) {
        ierr = MatDestroy(&D);
        CHKERRQ(ierr);
      }
      ierr = VecDestroy(&f);
      CHKERRQ(ierr);
      ierr = VecDestroy(&w);
      CHKERRQ(ierr);
      ierr = VecDestroy(&g);
      CHKERRQ(ierr);
      ierr = VecDestroy(&z);
      CHKERRQ(ierr);
    }
    end_time=clock();
    UpdateTime+=(end_time-start_time);

    /// ***************** END SOLVE AND UPDATE *******************
   } //end loop over vanka block
  } //loop over subdomain

  // *** Computational info ***
  cout << "Grid: "<<_msh->GetGridNumber()<< "      SOLVER TIME:              " <<
       static_cast<double>(SearchTime + AssemblyTime + SolveTime0 + SolveTime1 + SolveTime2 + UpdateTime)/ CLOCKS_PER_SEC<<
       "  ITS: " << its_A + its_C << endl;


//    cout <<  " SearchTime: " << static_cast<double>(SearchTime)/ CLOCKS_PER_SEC << endl;
//    cout << " AssemblyTime: " << static_cast<double>(AssemblyTime)/ CLOCKS_PER_SEC << endl;
//    cout << " SolverTime: " <<  static_cast<double>(SolveTime0 + SolveTime1 + SolveTime2)/ CLOCKS_PER_SEC << endl;
//    cout << " UpdateTime: " <<  static_cast<double>(UpdateTime)/ CLOCKS_PER_SEC << endl;

  return ierr;

}

//------------------------------------------------------------------------------------------------

int PetscLinearSolver::Vanka_Smoother(const vector <unsigned> &MGIndex, const vector <unsigned> &VankaIndex) {
  PetscErrorCode ierr;
  clock_t SearchTime=0, AssemblyTime=0, SolveTime0=0, UpdateTime=0;
  clock_t start_time, end_time;

  static PetscInt A_max=0;
  static PetscInt C_max=0;
  int its_0=0, its=0;
  unsigned int *end_ind = new unsigned int[MGIndex.size()];

  for (unsigned indexSol=0; indexSol<MGIndex.size(); indexSol++) {
    end_ind[indexSol]=_msh->_END_IND[_SolType[MGIndex[indexSol]]];
  }

  unsigned nvt=KKIndex[KKIndex.size()-1u];
  unsigned nel=_msh->GetElementNumber();

  if (indexa.size()<nvt || indexc.size()<nel) {
    indexa.resize(nvt);
    indexb.resize(nvt);
    indexc.resize(nel);
    indexd.resize(nel);
  }

  unsigned IndexaSize=50000;
  unsigned IndexbSize=50000;
  unsigned IndexcSize=10000;
  unsigned IndexdSize=10000;

  if (indexai.size()<IndexaSize) {
    indexai.resize(IndexaSize);
    indexbi.resize(IndexbSize);
    indexci.resize(IndexcSize);
    indexdi.resize(IndexdSize);
  }

  for (unsigned i=0; i<nvt; i++) {
    indexa[i]=IndexaSize;
    indexb[i]=IndexbSize;
  }
  for (unsigned i=0; i<nel; i++) {
    indexc[i]=IndexcSize;
    indexd[i]=IndexdSize;
  }

  // *** Start Vanka Block ***
  int grid=_msh->GetGridNumber();
  
 for (int gel=0; gel < static_cast <int> (nel); gel += _num_elem_vanka_block) {
    /// ***************** NODE/ELEMENT SERCH *******************
    start_time=clock();
    PetscInt Asize=0;
    PetscInt counterb=0;
    PetscInt Csize=0;
    PetscInt Dsize=0;
    for (int iel=gel; iel<gel+_num_elem_vanka_block && iel< static_cast <int> (nel); iel++) {
      for (unsigned i=0; i<_msh->el->GetElementDofNumber(iel,0); i++) {
        unsigned inode=_msh->el->GetElementVertexIndex(iel,i)-1u;
        unsigned nvei=_msh->el->GetVertexElementNumber(inode);
        const unsigned *pt_jel=_msh->el->GetVertexElementAddress(inode,0);
        for (unsigned j=0; j<nvei; j++) {
          unsigned jel=*(pt_jel++)-1u;
          //add elements for velocity to be solved
          if (indexc[jel]==IndexcSize) {
            indexci[Csize]=jel;
            indexc[jel]=Csize++;
            ///add node variables to be solved
            for (unsigned iind=0; iind<VankaIndex.size(); iind++) { 
              unsigned indexSol=VankaIndex[iind];
              const unsigned *pt_un=_msh->el->GetElementVertexAddress(jel,0);
              unsigned nvej=_msh->el->GetElementDofNumber(jel,end_ind[indexSol]);
              for (unsigned jj=0; jj<nvej; jj++) {
                unsigned jnode=*(pt_un++)-1u;
                if (indexa[jnode+KKIndex[indexSol]]==IndexaSize && 1.1<(*(*_Bdc)[MGIndex[indexSol]])(jnode)) {
                  indexai[Asize]=jnode+KKIndex[indexSol];
                  indexa[jnode+KKIndex[indexSol]]=Asize++;
                }
              }
            }
            for (unsigned jj=0; jj<_msh->el->GetElementDofNumber(jel,0); jj++) {
              unsigned jnode=_msh->el->GetElementVertexIndex(jel,jj)-1u;
              unsigned nvej=_msh->el->GetVertexElementNumber(jnode);
              const unsigned *pt_kel=_msh->el->GetVertexElementAddress(jnode,0);
              for (unsigned k=0; k<nvej; k++) {
                unsigned kel=*(pt_kel++)-1u;
                //add elements for all variables to be updated
                if (indexd[kel]==IndexdSize) {
                  indexdi[Dsize]=kel;
                  indexd[kel]=Dsize++;
                  ///add all variables to be updated
                  for (unsigned indexSol=0; indexSol<KKIndex.size()-1u; indexSol++) {
                    const unsigned *pt_un=_msh->el->GetElementVertexAddress(kel,0);
                    unsigned nvek=_msh->el->GetElementDofNumber(kel,end_ind[indexSol]);
                    for (unsigned kk=0; kk<nvek; kk++) {
                      unsigned knode=*(pt_un++)-1u;
                      if (indexb[knode+KKIndex[indexSol]]==IndexbSize && 0.1<(*(*_Bdc)[MGIndex[indexSol]])(knode)) {
                        indexbi[counterb]=knode+KKIndex[indexSol];
                        indexb[knode+KKIndex[indexSol]]=counterb++;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if (0==Asize) continue;

    PetscInt PBsize=Asize;
    for (PetscInt i=0; i<counterb; i++) {
      unsigned jnode=indexbi[i];
      if (indexa[jnode]==IndexaSize) {
        indexai[Asize]=jnode;
        indexa[jnode]=Asize++;
      }
      // *** reinitialize indexb
      indexb[jnode]=IndexbSize;
    }
    // *** re-initialize indeces
    for (PetscInt i=0; i<Asize; i++) {
      indexa[indexai[i]]=IndexaSize;
    }
    for (PetscInt i=0; i<Csize; i++) {
      indexc[indexci[i]]=IndexcSize;
    }
    for (PetscInt i=0; i<Dsize; i++) {
      indexd[indexdi[i]]=IndexdSize;
    }

    PetscInt PAmBsize=Asize-PBsize;

    //sort only for solved variables
    sort(indexai.begin(),indexai.begin()+PBsize);

    //sort only for updating variables but not solved
    sort(indexai.begin()+PBsize,indexai.begin()+Asize);

    if (Asize>A_max) A_max=Asize;
    if (Csize>C_max) C_max=Csize;

    end_time=clock();
    SearchTime+=(end_time-start_time);
    start_time=clock();
    /// ***************** END NODE/ELEMENT SERCH *******************

    /// ***************** ASSEMBLY *******************


    start_time=clock();
    Vec Pw, Pf, Pr;     /* exact solution, RHS*/
    Mat PA    ;         /* linear system matrix */

    // generate IS
    IS isA;
    ierr=ISCreateGeneral(PETSC_COMM_SELF,PBsize,&indexai[0],PETSC_USE_POINTER,&isA);
    CHKERRQ(ierr);

    // Solution Vector Pw
    ierr = VecCreateSeq(PETSC_COMM_SELF,PBsize,&Pw);
    CHKERRQ(ierr);
//     ierr = VecSetSizes(Pw,PETSC_DECIDE,PBsize);
//     CHKERRQ(ierr);
    ierr = VecSetFromOptions(Pw);
    CHKERRQ(ierr);

    // RHS Pf
    ierr = VecDuplicate(Pw,&Pf);
    CHKERRQ(ierr);

    PetscInt *indl=new PetscInt[PBsize];
    for (int i=0; i<PBsize; i++) indl[i]=i;
    PetscScalar *y=new PetscScalar [PBsize];
    const PetscInt *ind[1];

    ierr = ISGetIndices(isA,ind);
    CHKERRQ(ierr);
    ierr = VecGetValues(RES,PBsize,ind[0],y);
    ierr = VecSetValues(Pf,PBsize,indl,y,INSERT_VALUES);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(Pf);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(Pf);
    CHKERRQ(ierr);

    ierr = ISRestoreIndices(isA,ind);
    CHKERRQ(ierr);

    delete [] y;
    delete [] indl;


    // RHS Pr
    ierr = VecDuplicate(Pw,&Pr);
    CHKERRQ(ierr);

    //Matrix PA
    ierr = MatGetSubMatrix(KK,isA,isA,MAT_INITIAL_MATRIX,&PA);
    CHKERRQ(ierr);

    end_time=clock();
    AssemblyTime+=(end_time-start_time);

    /// ***************** END ASSEMBLY ******************

    /// ***************** SOLVE *******************
    start_time=clock();
    //set KSP and solve

    // Initialize the data structure with the matrix PA
    this->init(PA);

    // Solve the linear system
    ierr = KSPSolve(_ksp,Pf,Pw);
    CHKERRQ(ierr);

    // print information on ksp solver
    //ierr = KSPView(_ksp,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

    ierr = KSPGetIterationNumber(_ksp,&its);
    its_0 += its;

    end_time=clock();
    SolveTime0+=(end_time-start_time);
    /// ***************** END SOLVE ******************

    /// ***************** START UPDATE ******************

    start_time=clock();

    //update Residual and solution for the nodes that are solved for
    ierr = MatMult(PA,Pw,Pr);
    CHKERRQ(ierr);

    ierr = VecScale (Pr, -1.);
    CHKERRQ(ierr);
    PetscScalar *R[1];
    ierr = VecGetArray(Pr,R);
    CHKERRQ(ierr);

    PetscScalar *W[1];
    ierr = VecGetArray(Pw,W);
    CHKERRQ(ierr);

    ierr = ISGetIndices(isA,ind);
    CHKERRQ(ierr);

    ierr=VecSetValues(RES,PBsize,ind[0],R[0],ADD_VALUES);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(RES);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(RES);
    CHKERRQ(ierr);

    ierr=VecSetValues(EPS,PBsize,ind[0],W[0],ADD_VALUES);
    CHKERRQ(ierr);
    ierr = VecAssemblyBegin(EPS);
    CHKERRQ(ierr);
    ierr = VecAssemblyEnd(EPS);
    CHKERRQ(ierr);

    ierr = ISRestoreIndices(isA,ind);
    CHKERRQ(ierr);
    ierr = VecRestoreArray(Pr,R);
    CHKERRQ(ierr);

    ierr = VecRestoreArray(Pw,W);
    CHKERRQ(ierr);

    if (PAmBsize) { //update residual for the nodes that are not solved for
      IS isB;
      ierr=ISCreateGeneral(PETSC_COMM_SELF,PAmBsize,&indexai[PBsize],PETSC_USE_POINTER,&isB);
      CHKERRQ(ierr);

      Vec Ps;
      ierr = VecCreateSeq(PETSC_COMM_SELF,PAmBsize,&Ps);
      CHKERRQ(ierr);
//       ierr = VecSetSizes(Ps,PETSC_DECIDE,PAmBsize);
//       CHKERRQ(ierr);
      ierr = VecSetFromOptions(Ps);
      CHKERRQ(ierr);

      Mat PB;
      ierr = MatGetSubMatrix(KK,isB,isA,MAT_INITIAL_MATRIX,&PB);
      CHKERRQ(ierr);
      ierr = MatMult(PB,Pw,Ps);
      CHKERRQ(ierr);

      ierr = VecScale (Ps, -1.);
      CHKERRQ(ierr);
      PetscScalar *S[1];
      ierr = VecGetArray(Ps,S);
      CHKERRQ(ierr);
      ierr = ISGetIndices(isB,ind);
      CHKERRQ(ierr);

      ierr = VecSetValues(RES,PAmBsize,ind[0],S[0],ADD_VALUES);
      CHKERRQ(ierr);
      ierr = VecAssemblyBegin(RES);
      CHKERRQ(ierr);
      ierr = VecAssemblyEnd(RES);
      CHKERRQ(ierr);

      ierr = ISRestoreIndices(isB,ind);
      CHKERRQ(ierr);
      ierr = VecRestoreArray(Ps,S);
      CHKERRQ(ierr);
      ierr = VecDestroy(&Ps);
      CHKERRQ(ierr);
      ierr = MatDestroy(&PB);
      CHKERRQ(ierr);
      ierr = ISDestroy(&isB);
      CHKERRQ(ierr);
    }

    ierr = VecDestroy(&Pw);
    CHKERRQ(ierr);
    ierr = VecDestroy(&Pf);
    CHKERRQ(ierr);
    ierr = VecDestroy(&Pr);
    CHKERRQ(ierr);
    ierr = MatDestroy(&PA);
    CHKERRQ(ierr);
    ierr = ISDestroy(&isA);
    CHKERRQ(ierr);

    ierr = KSPDestroy(&_ksp);
    CHKERRQ(ierr);

    end_time=clock();
    UpdateTime+=(end_time-start_time);

    /// ***************** END SOLVE AND UPDATE *******************
  }
  
  // *** Computational info ***
  cout << "Grid: "<<_msh->GetGridNumber()<< "      SOLVER TIME:              " <<
       static_cast<double>(SearchTime + AssemblyTime + SolveTime0 + UpdateTime)/ CLOCKS_PER_SEC<<
       "  ITS: " << its_0 << endl;

  delete[] end_ind;

  return 1;
}


//-----------------------------------------------------------------------------------------------------


void PetscLinearSolver::clear() {
  if (this->initialized()) {
    this->_is_initialized = false;
    int ierr=0;
    ierr = KSPDestroy(&_ksp);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Mimic PETSc default solver and preconditioner
    this->_solver_type  = GMRES;
    if (!this->_preconditioner)      {
      int i;
      MPI_Comm_size(MPI_COMM_WORLD,&i);
      if (i == 1) this->_preconditioner_type = ILU_PRECOND;
      else this->_preconditioner_type = BLOCK_JACOBI_PRECOND;
    }
  }
}

// ==============================================================
void PetscLinearSolver::init() {
  // Initialize the data structures if not done so already.
  if (!this->initialized()) {
    this->_is_initialized = true;
    int ierr=0;
    // Create the linear solver context
    ierr = KSPCreate(MPI_COMM_WORLD, &_ksp);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Create the preconditioner context
    ierr = KSPGetPC(_ksp, &_pc);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Have the Krylov subspace method use our good initial guess rather than 0
    ierr = KSPSetInitialGuessNonzero(_ksp, PETSC_TRUE);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Set user-specified  solver and preconditioner types
    this->set_petsc_solver_type();
    // Set the options from user-input
    // Set runtime options, e.g., -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    //  These options will override those specified above as long as
    //  KSPSetFromOptions() is called _after_ any other customization  routines.
    ierr = KSPSetFromOptions(_ksp);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
    // NOT NECESSARY!!!!
    //ierr = PCSetFromOptions (_pc);
    //CHKERRABORT(MPI_COMM_WORLD,ierr);
// #endif

    // Notify PETSc of location to store residual history.
    // This needs to be called before any solves, since
    // it sets the residual history length to zero.  The default
    // behavior is for PETSc to allocate (internally) an array
    // of size 1000 to hold the residual norm history.
    ierr = KSPSetResidualHistory(_ksp,
                                 PETSC_NULL,   // pointer to the array which holds the history
                                 PETSC_DECIDE, // size of the array holding the history
                                 PETSC_TRUE);  // Whether or not to reset the history for each solve.
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);

    //If there is a preconditioner object we need to set the internal setup and apply routines
    if (this->_preconditioner) {
      PCShellSetContext(_pc,(void*)this->_preconditioner);
      PCShellSetSetUp(_pc,__libmesh_petsc_preconditioner_setup);
      PCShellSetApply(_pc,__libmesh_petsc_preconditioner_apply);
    }
  }
}


// ========================================================
void PetscLinearSolver::init(Mat& matrix, const bool pc_flag, const bool Schur) {
  // Initialize the data structures if not done so already.
  //I have to initialize many times the data structure therefore the next line must be commented
//   if (!this->initialized()) {
  this->_is_initialized = true;
  int ierr=0;
  int grid=_msh->GetGridNumber();

  // Since the matrix is sequential for us, we use as communicator PETSC_COMM_SELF (serial) 
  // instead of PETSC_COMM_WORLD (parallel)
  ierr = KSPCreate(PETSC_COMM_SELF,&_ksp);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = KSPSetOperators(_ksp,matrix,matrix,DIFFERENT_NONZERO_PATTERN);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = KSPSetInitialGuessNonzero(_ksp, PETSC_TRUE);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  if (!Schur || grid==0) {
    ierr = KSPSetTolerances(_ksp,_rtol,_abstol,_dtol,_maxits);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  } else {
    ierr = KSPSetTolerances(_ksp,_rtol_,_abstol_,_dtol_,_maxits_);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  ierr = KSPSetFromOptions(_ksp);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Set user-specified  solver and preconditioner types
  this->set_petsc_solver_type();

  if (!pc_flag) {
    ierr = KSPGetPC(_ksp,&_pc);
    CHKERRABORT(MPI_COMM_WORLD,ierr);  
    PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);
    PetscReal zero = 1.e-16;
    ierr = PCFactorSetZeroPivot(_pc,zero); CHKERRABORT(MPI_COMM_WORLD,ierr); 
    ierr = PCFactorSetShiftType(_pc,MAT_SHIFT_NONZERO); CHKERRABORT(MPI_COMM_WORLD,ierr); 
    
  }
}

// ========================================================
void PetscLinearSolver::init_schur(Mat& matrix) {
  // Initialize the data structures if not done so already.
  //I have to initialize many times the data structure therefore the next line must be commented
//   if (!this->initialized()) {
  this->_is_initialized = true;
  int ierr=0;

  ierr = MatSchurComplementGetKSP(matrix,&_ksp2);
  ierr = KSPSetInitialGuessNonzero(_ksp2, PETSC_TRUE);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  ierr = KSPSetTolerances(_ksp2,_rtol,_abstol,_dtol,_maxits);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = KSPSetFromOptions(_ksp2);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  // Set user-specified  solver and preconditioner types
  this->set_petsc_solver_type2();

  ierr = KSPGetPC(_ksp2, &_pc);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);
  PetscReal zero = 1.e-16;
  PCFactorSetZeroPivot(_pc,zero);
  PCFactorSetShiftType(_pc,MAT_SHIFT_NONZERO);

//   }
}

// ========================================================
void PetscLinearSolver::init(Mat& matrix) {
  // Initialize the data structures if not done so already.
  if (!this->initialized())    {
    this->_is_initialized = true;
    int ierr=0;
    // Create the linear solver context
    ierr = KSPCreate(MPI_COMM_WORLD, &_ksp);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    //ierr = PCCreate (MPI_COMM_WORLD, &_pc); CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Create the preconditioner context
    ierr = KSPGetPC(_ksp, &_pc);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Set operators. The input matrix works as the preconditioning matrix
    ierr = KSPSetOperators(_ksp, matrix, matrix, SAME_NONZERO_PATTERN);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Have the Krylov subspace method use our good initial guess rather than 0
    ierr = KSPSetInitialGuessNonzero(_ksp, PETSC_TRUE);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Set user-specified  solver and preconditioner types
    this->set_petsc_solver_type();
    
    if (!this->same_preconditioner)  {
      ierr = KSPSetOperators(_ksp, KK, KK, SAME_NONZERO_PATTERN);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    } else {
      ierr = KSPSetOperators(_ksp, KK, KK, SAME_PRECONDITIONER);
      CHKERRABORT(MPI_COMM_WORLD,ierr);
    }

    // Set the tolerances for the iterative solver.  Use the user-supplied
    // tolerance for the relative residual & leave the others at default values.
    ierr = KSPSetTolerances(_ksp,_rtol,_abstol,_dtol,_maxits);
    //   ierr = KSPSetTolerances(_ksp, , PETSC_DEFAULT,PETSC_DEFAULT, max_its);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    
    // Set the options from user-input
    // Set runtime options, e.g., -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    //  These options will override those specified above as long as
    //  KSPSetFromOptions() is called _after_ any other customization  routines.
    ierr = KSPSetFromOptions(_ksp);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    // Not sure if this is necessary, or if it is already handled by KSPSetFromOptions?
    //ierr = PCSetFromOptions (_pc);CHKERRABORT(MPI_COMM_WORLD,ierr);

    // Notify PETSc of location to store residual history.
    // This needs to be called before any solves, since
    // it sets the residual history length to zero.  The default
    // behavior is for PETSc to allocate (internally) an array
    // of size 1000 to hold the residual norm history.
    ierr = KSPSetResidualHistory(_ksp,
                                 PETSC_NULL,   // pointer to the array which holds the history
                                 PETSC_DECIDE, // size of the array holding the history
                                 PETSC_TRUE);  // Whether or not to reset the history for each solve.
    CHKERRABORT(MPI_COMM_WORLD,ierr);

    PetscPreconditioner::set_petsc_preconditioner_type(this->_preconditioner_type,_pc);
    PetscReal zero = 1.e-16;
    PCFactorSetZeroPivot(_pc,zero);
    PCFactorSetShiftType(_pc,MAT_SHIFT_NONZERO);
//     if (this->_preconditioner) {
//       this->_preconditioner->set_matrix(*matrix);
//       PCShellSetContext(_pc,(void*)this->_preconditioner);
//       PCShellSetSetUp(_pc,__libmesh_petsc_preconditioner_setup);
//       PCShellSetApply(_pc,__libmesh_petsc_preconditioner_apply);
//     }
  }
}

// ========================================================
std::pair< int, double> PetscLinearSolver::solve() {

  

  int ierr=0;
  int its=0;
  int max_its =100; //= static_cast<int>(m_its);
  PetscReal final_resid=0.;
  
  //ierr = MatAssemblyBegin(KK,MAT_FINAL_ASSEMBLY);CHKERRABORT(MPI_COMM_WORLD,ierr);
  //ierr = MatAssemblyEnd(KK,MAT_FINAL_ASSEMBLY);CHKERRABORT(MPI_COMM_WORLD,ierr);
 //  MatZeroRows(KK,DrchKKdofs.size(),&DrchKKdofs[0],1.e20,PETSC_NULL,PETSC_NULL);

//    cout << endl << endl << "after BDC " << endl;
//           PetscViewer viewer;
//         PetscViewerSetFormat(viewer,PETSC_VIEWER_DEFAULT);
// // //       PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,700,700,&viewer);
//         MatView(KK,viewer);
  
  vector <PetscScalar>  value(DrchKKdofs.size());
  for(int i=0;i<DrchKKdofs.size();i++){
    PetscInt row=DrchKKdofs[i];
    ierr = MatGetValues(KK,1,&row,1,&row,&value[i]);
  }
   
  PetscScalar penalty=1.0e20;
  for(int i=0;i<DrchKKdofs.size();i++){
    PetscInt row=DrchKKdofs[i];
    ierr = MatSetValuesBlocked(KK,1,&row,1,&row,&penalty,INSERT_VALUES);  CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = MatAssemblyBegin(KK,MAT_FINAL_ASSEMBLY);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = MatAssemblyEnd(KK,MAT_FINAL_ASSEMBLY);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  
//             PetscViewer viewer;
//         PetscViewerSetFormat(viewer,PETSC_VIEWER_DEFAULT);
// // //       PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL,PETSC_NULL,0,0,700,700,&viewer);
//         MatView(KK,viewer);
  
  
  
//   vector< PetscScalar > zero(DrchKKdofs.size(),0.);
//   VecSetValues(RES,DrchKKdofs.size(),&DrchKKdofs[0],&zero[0],INSERT_VALUES);
//   
//   ierr = VecAssemblyBegin(RES);CHKERRABORT(MPI_COMM_WORLD,ierr);
//   ierr = VecAssemblyEnd(RES);CHKERRABORT(MPI_COMM_WORLD,ierr);
//   
//    if(_nprocs!=1) {
//     ierr = VecGhostUpdateBegin(RES,INSERT_VALUES,SCATTER_FORWARD);CHKERRABORT(MPI_COMM_WORLD,ierr);
//     ierr = VecGhostUpdateEnd(RES,INSERT_VALUES,SCATTER_FORWARD);CHKERRABORT(MPI_COMM_WORLD,ierr);
//   }
  
  
//   PetscReal value;
//   VecNorm(RES,NORM_2,&value);
//   cout << "res before solver : " << value << endl;
//   
//   VecNorm(EPS,NORM_2,&value);
//   cout << "eps before solver : " << value << endl;
  
//   // Close the matrices and vectors in case this wasn't already done.
//   matrix->close();
//   precond->close();
//   solution->close();
//   rhs->close();
// //   // If matrix != precond, then this means we have specified a
// //   // special preconditioner, so reset preconditioner type to PCMAT.
// //   if (matrix != precond)
// //     {
// //       this->_preconditioner_type = USER_PRECOND;
// //       this->set_petsc_preconditioner_type ();
// //     }
// 
//   if (this->_preconditioner) {
//     this->_preconditioner->set_matrix(matrix_in);
//   }
// 
//   // 2.2.1 & newer style
//   // Set operators. The input matrix works as the preconditioning matrix
// 
// 
//   if ( this->_preconditioner_type==MCC_PRECOND || this->_preconditioner_type==ICC_PRECOND ) {
// //       /* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
//     ierr = MatSetOption(matrix->mat(),MAT_SYMMETRIC,PETSC_TRUE);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     ierr = MatSetOption(matrix->mat(),MAT_SPD,PETSC_TRUE);
//     CHKERRABORT(MPI_COMM_WORLD,ierr); /* set MUMPS id%SYM=1 */
//   }
// 
//   if ( this->_preconditioner_type==MLU_PRECOND ) {
//     PetscInt ival,icntl;
// 
//     ierr = KSPSetOperators(_ksp, KK, KK, DIFFERENT_NONZERO_PATTERN);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//     ierr = PCFactorSetUpMatSolverPackage(_pc);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);   /* call MatGetFactor() to create F */
//     ierr = PCFactorGetMatrix(_pc,&KK);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     icntl=7;
//     ival = 2;
//     ierr = MatMumpsSetIcntl(KK,icntl,ival);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//   }
// 
//   else if ( this->_preconditioner_type==MCC_PRECOND ) {
//     PetscInt ival,icntl;
// 
//     ierr = KSPSetOperators(_ksp, matrix->mat(), precond_, DIFFERENT_NONZERO_PATTERN);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//     ierr = PCFactorSetUpMatSolverPackage(_pc);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);   /* call MatGetFactor() to create F */
//     ierr = PCFactorGetMatrix(_pc,&precond_);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     icntl=7;
//     ival = 2;
//     ierr = MatMumpsSetIcntl(precond_,icntl,ival);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//   else {
// 
  
  this->init(KK);  //as before
  
  VecZeroEntries(EPSC);
  
  // Solve the linear system
  ierr = KSPSolve(_ksp, RES, EPSC);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  
  ierr = VecAXPY(EPS,1.,EPSC);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  
  if(_msh->_nprocs!=1) {
    ierr = VecGhostUpdateBegin(EPS,INSERT_VALUES,SCATTER_FORWARD);CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostUpdateEnd(EPS,INSERT_VALUES,SCATTER_FORWARD);CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  
  ierr = VecScale(EPSC,-1.);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  
  for(int i=0;i<DrchKKdofs.size();i++){
    PetscInt row=DrchKKdofs[i];
    ierr = MatSetValuesBlocked(KK,1,&row,1,&row,&value[i],INSERT_VALUES);  CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = MatAssemblyBegin(KK,MAT_FINAL_ASSEMBLY);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = MatAssemblyEnd(KK,MAT_FINAL_ASSEMBLY);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  
  
   
  
  /*
  ierr=MatMult(KK,EPSC,RESC);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  
  vector< PetscScalar > zero(DrchKKdofs.size(),0.);
  VecSetValues(RESC,DrchKKdofs.size(),&DrchKKdofs[0],&zero[0],INSERT_VALUES);
  ierr = VecAssemblyBegin(RESC);CHKERRABORT(MPI_COMM_WORLD,ierr);
  ierr = VecAssemblyEnd(RESC);CHKERRABORT(MPI_COMM_WORLD,ierr);
  
  ierr = VecAXPY(RES,1.,RESC);*/
  
  
  ierr = MatMultAdd(KK,EPSC,RES,RES);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  
  if(_msh->_nprocs!=1) {
    ierr = VecGhostUpdateBegin(RES,INSERT_VALUES,SCATTER_FORWARD);CHKERRABORT(MPI_COMM_WORLD,ierr);
    ierr = VecGhostUpdateEnd(RES,INSERT_VALUES,SCATTER_FORWARD);CHKERRABORT(MPI_COMM_WORLD,ierr);
  }
  
//  VecZeroEntries(EPSC);
 
//   PetscReal value;
//   VecNorm(RES,NORM_2,&value);
//   cout << "res after solver : " << value << endl;
//   
//   VecNorm(EPS,NORM_2,&value);
//   cout << "eps after solver : " << value << endl;
 
  // Get the number of iterations required for convergence
  ierr = KSPGetIterationNumber(_ksp, &its);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  
  // Get the norm of the final residual to return to the user.
  ierr = KSPGetResidualNorm(_ksp, &final_resid);
  CHKERRABORT(MPI_COMM_WORLD,ierr);

  // STOP_LOG("solve()", "PetscLinearSolver");
  return std::make_pair(its, final_resid);
  
}


// // ========================================================
// std::pair< int, double> PetscLinearSolver::solve(SparseMatrix&  matrix_in,
//     SparseMatrix&  precond_in,  NumericVector& solution_in,  NumericVector& rhs_in,
//     const double tol,   const  int m_its) {
// 
// //   START_LOG("solve()", "PetscLinearSolver");
//   // Make sure the data passed in are really of Petsc types
//   PetscMatrix* matrix   = libmeshM_cast_ptr<PetscMatrix*>(&matrix_in);
//   PetscMatrix* precond  = libmeshM_cast_ptr<PetscMatrix*>(&precond_in);
//   PetscVector* solution = libmeshM_cast_ptr<PetscVector*>(&solution_in);
//   PetscVector* rhs      = libmeshM_cast_ptr<PetscVector*>(&rhs_in);
// 
//   Mat precond_ = precond->mat();
// 
//   this->init(matrix);  //as before
// 
// //   this->init();
// 
//   int ierr=0;
//   int its=0, max_its = static_cast<int>(m_its);
//   PetscReal final_resid=0.;
//   // Close the matrices and vectors in case this wasn't already done.
//   matrix->close();
//   precond->close();
//   solution->close();
//   rhs->close();
// //   // If matrix != precond, then this means we have specified a
// //   // special preconditioner, so reset preconditioner type to PCMAT.
// //   if (matrix != precond)
// //     {
// //       this->_preconditioner_type = USER_PRECOND;
// //       this->set_petsc_preconditioner_type ();
// //     }
// 
//   if (this->_preconditioner) {
//     this->_preconditioner->set_matrix(matrix_in);
//   }
// 
//   // 2.2.1 & newer style
//   // Set operators. The input matrix works as the preconditioning matrix
// 
// 
//   if ( this->_preconditioner_type==MCC_PRECOND || this->_preconditioner_type==ICC_PRECOND ) {
// //       /* A is symmetric. Set symmetric flag to enable ICC/Cholesky preconditioner */
//     ierr = MatSetOption(matrix->mat(),MAT_SYMMETRIC,PETSC_TRUE);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     ierr = MatSetOption(matrix->mat(),MAT_SPD,PETSC_TRUE);
//     CHKERRABORT(MPI_COMM_WORLD,ierr); /* set MUMPS id%SYM=1 */
//   }
// 
//   if ( this->_preconditioner_type==MLU_PRECOND ) {
//     PetscInt ival,icntl;
// 
//     ierr = KSPSetOperators(_ksp, matrix->mat(), precond_, DIFFERENT_NONZERO_PATTERN);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//     ierr = PCFactorSetUpMatSolverPackage(_pc);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);   /* call MatGetFactor() to create F */
//     ierr = PCFactorGetMatrix(_pc,&precond_);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     icntl=7;
//     ival = 2;
//     ierr = MatMumpsSetIcntl(precond_,icntl,ival);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//   }
// 
//   else if ( this->_preconditioner_type==MCC_PRECOND ) {
//     PetscInt ival,icntl;
// 
//     ierr = KSPSetOperators(_ksp, matrix->mat(), precond_, DIFFERENT_NONZERO_PATTERN);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//     ierr = PCFactorSetUpMatSolverPackage(_pc);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);   /* call MatGetFactor() to create F */
//     ierr = PCFactorGetMatrix(_pc,&precond_);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
//     icntl=7;
//     ival = 2;
//     ierr = MatMumpsSetIcntl(precond_,icntl,ival);
//     CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//   } else {
// 
//     if (!this->same_preconditioner)  {
//       ierr = KSPSetOperators(_ksp, matrix->mat(), precond_, SAME_NONZERO_PATTERN);
//       CHKERRABORT(MPI_COMM_WORLD,ierr);
//     } else {
//       ierr = KSPSetOperators(_ksp, matrix->mat(), precond_, SAME_PRECONDITIONER);
//       CHKERRABORT(MPI_COMM_WORLD,ierr);
//     }
//   }
// 
//   // Set the tolerances for the iterative solver.  Use the user-supplied
//   // tolerance for the relative residual & leave the others at default values.
//   ierr = KSPSetTolerances(_ksp, tol, PETSC_DEFAULT,PETSC_DEFAULT, max_its);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
// 
//   // Solve the linear system
// 
// // PetscLogEvent USER_EVENT;
// // PetscLogDouble user_event_flops;
// // PetscLogEventRegister("User event",0,&USER_EVENT);
// // PetscLogEventBegin(USER_EVENT,0,0,0,0);
// 
// 
//   ierr = KSPSolve(_ksp, rhs->vec(), solution->vec());
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//   // PetscLogFlops(user_event_flops);
//   // PetscLogEventEnd(USER_EVENT,0,0,0,0);
// 
//   // Get the number of iterations required for convergence
//   ierr = KSPGetIterationNumber(_ksp, &its);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
//   // Get the norm of the final residual to return to the user.
//   ierr = KSPGetResidualNorm(_ksp, &final_resid);
//   CHKERRABORT(MPI_COMM_WORLD,ierr);
// 
//   // STOP_LOG("solve()", "PetscLinearSolver");
//   return std::make_pair(its, final_resid);
// }


// =========================================================================
void PetscLinearSolver::get_residual_history(std::vector<double>& hist) {
  int ierr = 0;
  int its  = 0;
  // Fill the residual history vector with the residual norms
  // Note that GetResidualHistory() does not copy any values, it
  // simply sets the pointer p.  Note that for some Krylov subspace
  // methods, the number of residuals returned in the history
  // vector may be different from what you are expecting.  For
  // example, TFQMR returns two residual values per iteration step.
  PetscReal* p;
  ierr = KSPGetResidualHistory(_ksp, &p, &its);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Check for early return
  if (its == 0) return;
  // Create space to store the result
  hist.resize(its);
  // Copy history into the vector provided by the user.
  for (int i=0; i<its; ++i) {
    hist[i] = *p;
    p++;
  }
}

// ======================================================
double PetscLinearSolver::get_initial_residual() {
  int ierr = 0;
  int its  = 0;
  // Fill the residual history vector with the residual norms
  // Note that GetResidualHistory() does not copy any values, it
  // simply sets the pointer p.  Note that for some Krylov subspace
  // methods, the number of residuals returned in the history
  // vector may be different from what you are expecting.  For
  // example, TFQMR returns two residual values per iteration step.
  PetscReal* p;
  ierr = KSPGetResidualHistory(_ksp, &p, &its);
  CHKERRABORT(MPI_COMM_WORLD,ierr);
  // Check no residual history
  if (its == 0) {
    std::cerr << "No iterations have been performed, returning 0." << std::endl;
    return 0.;
  }
  // Otherwise, return the value pointed to by p.
  return *p;
}

// =================================================
void PetscLinearSolver::set_petsc_solver_type() {
  int ierr = 0;
  switch (this->_solver_type) {
  case CG:
    ierr = KSPSetType(_ksp, (char*) KSPCG);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case CR:
    ierr = KSPSetType(_ksp, (char*) KSPCR);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case CGS:
    ierr = KSPSetType(_ksp, (char*) KSPCGS);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case BICG:
    ierr = KSPSetType(_ksp, (char*) KSPBICG);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case TCQMR:
    ierr = KSPSetType(_ksp, (char*) KSPTCQMR);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case TFQMR:
    ierr = KSPSetType(_ksp, (char*) KSPTFQMR);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case LSQR:
    ierr = KSPSetType(_ksp, (char*) KSPLSQR);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case BICGSTAB:
    ierr = KSPSetType(_ksp, (char*) KSPBCGS);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case MINRES:
    ierr = KSPSetType(_ksp, (char*) KSPMINRES);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case GMRES:
    ierr = KSPSetType(_ksp, (char*) KSPGMRES);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case RICHARDSON:
    ierr = KSPSetType(_ksp, (char*) KSPRICHARDSON);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case CHEBYSHEV:
    ierr = KSPSetType(_ksp, (char*) KSPCHEBYSHEV);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case PREONLY:
    ierr = KSPSetType(_ksp, (char*) KSPPREONLY);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  default:
    std::cerr << "ERROR:  Unsupported PETSC Solver: "
              << this->_solver_type               << std::endl
              << "Continuing with PETSC defaults" << std::endl;
  }
}


// =================================================
void PetscLinearSolver::set_petsc_solver_type2() {
  int ierr = 0;
  switch (this->_solver_type) {
  case CG:
    ierr = KSPSetType(_ksp2, (char*) KSPCG);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case CR:
    ierr = KSPSetType(_ksp2, (char*) KSPCR);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case CGS:
    ierr = KSPSetType(_ksp2, (char*) KSPCGS);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case BICG:
    ierr = KSPSetType(_ksp2, (char*) KSPBICG);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case TCQMR:
    ierr = KSPSetType(_ksp2, (char*) KSPTCQMR);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case TFQMR:
    ierr = KSPSetType(_ksp2, (char*) KSPTFQMR);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case LSQR:
    ierr = KSPSetType(_ksp2, (char*) KSPLSQR);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case BICGSTAB:
    ierr = KSPSetType(_ksp2, (char*) KSPBCGS);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case MINRES:
    ierr = KSPSetType(_ksp2, (char*) KSPMINRES);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case GMRES:
    ierr = KSPSetType(_ksp2, (char*) KSPGMRES);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case RICHARDSON:
    ierr = KSPSetType(_ksp2, (char*) KSPRICHARDSON);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case CHEBYSHEV:
    ierr = KSPSetType(_ksp2, (char*) KSPCHEBYSHEV);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  case PREONLY:
    ierr = KSPSetType(_ksp2, (char*) KSPPREONLY);
    CHKERRABORT(MPI_COMM_WORLD,ierr);
    return;
  default:
    std::cerr << "ERROR:  Unsupported PETSC Solver: "
              << this->_solver_type               << std::endl
              << "Continuing with PETSC defaults" << std::endl;
  }
}

// =======================================================
void PetscLinearSolver::print_converged_reason() {
  KSPConvergedReason reason;
  KSPGetConvergedReason(_ksp, &reason);
  std::cout << "Linear solver convergence/divergence reason: " << KSPConvergedReasons[reason] << std::endl;
//   //  KSP_CONVERGED_RTOL (residual 2-norm decreased by a factor of rtol, from 2-norm of right hand side)
//   //  KSP_CONVERGED_ATOL (residual 2-norm less than abstol)
//   //  KSP_CONVERGED_ITS (used by the preonly preconditioner that always uses ONE iteration)
//   //  KSP_CONVERGED_STEP_LENGTH
//   //  KSP_DIVERGED_ITS  (required more than its to reach convergence)
//   //  KSP_DIVERGED_DTOL (residual norm increased by a factor of divtol)
//   //  KSP_DIVERGED_NAN (residual norm became Not-a-number likely do to 0/0)
//   //  KSP_DIVERGED_BREAKDOWN (generic breakdown in method)
//   switch (reason) {
//   case KSP_CONVERGED_RTOL: {
//     std::cout << "Linear solver converged, relative tolerance reached." << std::endl;
//     break;
//   }
//   case KSP_CONVERGED_ATOL: {
//     std::cout << "Linear solver converged, absolute tolerance reached." << std::endl;
//     break;
//   }
//   // Divergence
//   case KSP_DIVERGED_ITS: {
//     std::cout << "Linear solver diverged, max no. of iterations reached." << std::endl;
//     break;
//   }
//   case KSP_DIVERGED_DTOL: {
//     std::cout << "Linear solver diverged, residual norm increase by dtol (default 1.e5)." << std::endl;
//     break;
//   }
// //   case KSP_DIVERGED_NAN: {
// //     std::cout << "Linear solver diverged, residual norm is NaN." << std::endl;
// //     break;
// //   }
//   case KSP_DIVERGED_BREAKDOWN: {
//     std::cout << "Linear solver diverged, generic breakdown in the method." << std::endl;
//     break;
//   }
//   default: {
//     std::cout << "Unknown/unsupported con(di)vergence reason: " << reason << std::endl;
//   }
//   }
}

#endif // #ifdef LIBMESH_HAVE_PETSC