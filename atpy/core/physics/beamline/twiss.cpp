
#ifndef _TWISS_CPP_
#define _TWISS_CPP_

#include "twiss.h"


Matrix2d _symplectic_conjugate_2by2(Matrix2d& m){
    Matrix2d tmp;
    tmp<<m(1,1), -m(0,1), -m(1,0), m(0,0);
    return tmp;
}

int calculate_coupled_period_twiss(double* mat, double* tws){
    Map<Matrix6d > M(mat,6,6 );
    double t1, t2, Delta, _sign;
    
    Matrix2d A, B, C, D, Bbar_and_C, R, mat_X,mat_Y, tmp;
    // Matrix<double,4,4,RowMajor> eye44;
    // Matrix<double,4,1> eta;

    A.block<2,2>(0,0)=M.block<2,2>(0,0);
    B.block<2,2>(0,0)=M.block<2,2>(0,2);
    C.block<2,2>(0,0)=M.block<2,2>(2,0);
    D.block<2,2>(0,0)=M.block<2,2>(2,2);

    Bbar_and_C = _symplectic_conjugate_2by2(B) + C;
    t1 = 0.5*(A.trace() - D.trace()) ;
    Delta = t1*t1 + Bbar_and_C.determinant() ;
    if(Delta < 0.0){
        return -1;
    }
    _sign = t1>0 ? -1. : 1. ;
    t2=fabs(t1)+sqrt(Delta );
    if(fabs(t2)<1e-14 ){
        R.fill(0.) ;
    }
    else{
        R=Bbar_and_C*(_sign/t2);
    }
    
	mat_X = A-B*R;
	mat_Y = D+C*_symplectic_conjugate_2by2(R);

    if(mat_X.determinant()<0.99 || mat_Y.determinant()<0.99  ){
        return -2;
    }

    double cmux, cmuy, smux,smuy;
    double betax,betay,alphax,alphay,gammax,gammay;
	cmux=0.5*(mat_X(1,1)+mat_X(0,0) );
	cmuy=0.5*(mat_Y(1,1)+mat_Y(0,0) );
    if(fabs(cmux)>1 || fabs(cmuy)>1 ){
        if(fabs(cmux)>1){
            // let unstable violation at the unit of 1e2, so that the valid tune info can be obtained with floor.
            tws[NUX] = 1e6+1e2*floor((fabs(cmux)-1)/0.01) ;
        }
        else{
            smux=sqrt(1.0 -cmux*cmux)*SIGN(mat_X(0,1) );
            tws[NUX]=(smux>0? acos(cmux): PIx2-acos(cmux) )/PIx2;
        }
        if(fabs(cmuy)>1){
            // let unstable violation at the unit of 1e2, so that the valid tune info can be obtained with floor.
            tws[NUY] = 1e6+1e2*floor((fabs(cmuy)-1)/0.01) ;
        }
        else{
            smuy=sqrt(1.0 -cmuy*cmuy)*SIGN(mat_Y(0,1) );
            tws[NUY]=(smuy>0? acos(cmuy): PIx2-acos(cmuy) )/PIx2;
        }
        return -3;
    }
    else{
        smux=sqrt(1.0 -cmux*cmux)*SIGN(mat_X(0,1) );
        smuy=sqrt(1.0 -cmuy*cmuy)*SIGN(mat_Y(0,1) );
    }

	betax=mat_X(0,1)/smux;
	gammax=-mat_X(1,0)/smux;
	betay=mat_Y(0,1)/smuy;
	gammay=-mat_Y(1,0)/smuy;

	alphax=0.5*(mat_X(0,0)-mat_X(1,1) )/smux ;
	alphay=0.5*(mat_Y(0,0)-mat_Y(1,1) )/smuy ;
    // eye44.setIdentity(4,4);
	Map<Matrix<double,4,1>>(&tws[ETAX],4,1)=( Matrix4d::Identity(4,4) -M.block<4,4>(0,0) ).inverse()*M.block<4,1>(0,5) ;
    // COORD      ,BETAX       ,ALPHAX     ,GAMMAX     ,BETAY      ,
    //     ALPHAY     ,GAMMAY      ,ETAX       ,ETAPX      ,NUX        ,NUY
    tws[BETAX]=betax;
    tws[ALPHAX]=alphax;
    tws[GAMMAX]=gammax;
    tws[BETAY]=betay;
    tws[ALPHAY]=alphay;
    tws[GAMMAY]=gammay;
    // tws[ETAX] = eta[0]; tws[ETAPX] = eta[1]; tws[ETAY] = eta[2]; tws[ETAPY] = eta[3];
    tws[NUX]=(mat_X(0,1)>0? acos(cmux): PIx2-acos(cmux) )/PIx2;
    tws[NUY]=(mat_Y(0,1)>0? acos(cmuy): PIx2-acos(cmuy) )/PIx2;
    tws[TWSMODE]=1;
    return 0;
}


int propagate_twiss(double* tout, double* tin, double* mat){


    Map<Matrix6d > M(mat,6,6 );
    double  t;
    
    // std::cout<<"M: "<<endl<<M<<endl;

    Matrix2d A, B, C, D, R, mat_X, mat_Y, R0,_R0;

    A.block<2,2>(0,0)=M.block<2,2>(0,0);
    B.block<2,2>(0,0)=M.block<2,2>(0,2);
    C.block<2,2>(0,0)=M.block<2,2>(2,0);
    D.block<2,2>(0,0)=M.block<2,2>(2,2);

	R0=Map<Matrix2d >(&tin[R1], 2,2 );
	_R0=_symplectic_conjugate_2by2(R0);
    double mode;
	if(tin[TWSMODE] == 1.0){
		mat_X=A-B*R0;
		t=mat_X.determinant();
        if(t>0.1){
            R=(D*R0-C)*_symplectic_conjugate_2by2(mat_X);
            R/=t;
            mat_X/=sqrt(t);
            mat_Y=D+C*_R0;
            mat_Y/=sqrt(mat_Y.determinant() );
            mode=1.0;
        }
        else{
            mat_X=C-D*R0;
            mat_X/=sqrt(mat_X.determinant() );
            mat_Y=B+A*_R0;
            t=mat_Y.determinant() ;
            R=-(D+C*_R0)*_symplectic_conjugate_2by2(mat_Y);
            R/=t;
            mat_Y/=sqrt(t);
            mode=2.0;
        }
    }
	else if (tin[TWSMODE] == 2.0){
		mat_X=B+A*_R0;
		t=mat_X.determinant();
        if( t>0.1){
            R=-(D+C*_R0)*_symplectic_conjugate_2by2(mat_X) ;
            R/=t ;
            mat_X/=sqrt(t) ;
            mat_Y=C-D*R0 ;
            mat_Y/=sqrt(mat_Y.determinant()) ;
            mode=1.0 ;
        }
        else{
            mat_X=D+C*_R0 ;
            mat_X/=sqrt(mat_X.determinant() ) ;
            mat_Y=A-B*R0 ;
            t=mat_Y.determinant() ;
            R=(D*R0-C)*_symplectic_conjugate_2by2(mat_Y) ;
            R/=t ;
            mat_Y/=sqrt(t) ;
            mode=2.0 ;
        }
    }
	else{
		// throw(AssertionError("Mode should be integer 1 or 2."))
		// println(stderr,"Invalid mode.")
		return -1;
    }
	tout[R1]=R(0,0);
	tout[R2]=R(0,1);
	tout[R3]=R(1,0);
	tout[R4]=R(1,1);
    tout[NUX] = _propagate_decoupled_twiss_2x2(&tout[BETAX],&tin[BETAX], mat_X);
    tout[NUY] = _propagate_decoupled_twiss_2x2(&tout[BETAY],&tin[BETAY], mat_Y);
    // std::cout<<mat_X<<endl<<mat_Y<<endl;
    Map<Matrix<double,4,1>>(&tout[ETAX],4,1) = M.block<4,4>(0,0)*Map<Matrix<double,4,1>>(&tin[ETAX],4,1) + M.block<4,1>(0,5);
    tout[TWSMODE]=mode;
    return 0;
}

double  _propagate_decoupled_twiss_2x2(double* tout,double* tin, Matrix2d& mat){
    double m11,m12,m21,m22;
    double smux,cmux;
    m11=mat(0,0);
    m12=mat(0,1);
    m21=mat(1,0);
    m22=mat(1,1);
    tout[0] = m11*m11*tin[0]  - 2*m11*m12*tin[1] + m12*m12*tin[2];
    tout[1] = -m11*m21*tin[0] +  (1.0+2*m12*m21)*tin[1] -m12*m22*tin[2];
    tout[2] = m21*m21*tin[0]  - 2*m21*m22*tin[1] + m22*m22*tin[2];
    smux=m12/sqrt(tout[0]*tin[0]) ;
	cmux=m11*sqrt(tin[0]/tout[0])-tin[1]*smux ;

    // std::cout<<smux<<", "<<cmux<<endl;

    return (smux<0?  PIx2-acos(cmux):acos(cmux) )/PIx2; //consider smux==0
    // return fmod(atan2(smux,cmux), PIx2 );
}

#endif
