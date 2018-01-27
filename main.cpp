#include <math.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "TBody.h"

using namespace std;


const double AU=1.495978706910000E+08;
const double SunGm=2.959122082855911025e-04;


double _mu=2.959122082855911025e-04;
double _InitialElem[6]={1.5,0.01,1,1,1,1};
double _TC=0.00000000000000000000000000;

TBody _TBD;

TVector CalcuAcc(TVector R)
{ TVector A;
  double rr=R.Length();

  A=R*(-_mu/(rr*rr*rr));

  return A;
}



double _FT_back[3][15];

void CopyFT(int endk)
{ for(int j=0;j<3;j++) for(int k=-endk;k<=endk;k++) _FT_back[j][k+7]=_TBD.f[j][k+7];
}

double CmpFT(int endk)
{ double sum=0.00000000000000000000000000000;
  for(int j=0;j<3;j++) for(int k=-endk;k<=endk;k++)
    { if(_FT_back[j][k+7]==0.00000000000000000000000)
         sum+=fabs(_TBD.f[j][k+7]);
      else
         sum+=fabs((_TBD.f[j][k+7]-_FT_back[j][k+7])/_FT_back[j][k+7]);
    }
  return sum;
}

void Initia(int order,double h)
{ TVector R0,V0,A,R,V,A0,R1,V1,A1,Rf1,Vf1,Af1,Af1back;
  double h2,epslong;
  int ord;

  _TBD.h=h;
  h2=h;
  epslong=1.0e-25;
  _TC=0.0000000000000000000000;

  ElemToCor(_mu,0.0,_InitialElem,R0,V0);
  A0=CalcuAcc(R0);
  R=R0+V0*h+A0*(h*h/2); A1=CalcuAcc(R);
  R=R0*2-R;             Af1=CalcuAcc(R);

  _TBD.WriteF(-1,Af1);   _TBD.WriteF(0,A0);    _TBD.WriteF(1,A1);


  for(int ord=2;ord<=12;ord+=2)
    { do{ CopyFT(ord/2);

          _TBD.CalcuCentralSf3(ord,R0,V0);
          for(int k=1;k<=ord/2;k++)
            { _TBD.CalcuIntXv(k,ord,k,R,V);   A=CalcuAcc(R);
              _TBD.WriteF(k,A);
              _TBD.CalcuCentralSf3(ord,R0,V0);
              _TBD.CalcuSfAccCentral3(ord/2+1);
            }

          for(int k=-1;k>=-ord/2;k--)
            { _TBD.CalcuIntXv(k,ord,k,R,V);   A=CalcuAcc(R);
              _TBD.WriteF(k,A);
              _TBD.CalcuCentralSf3(ord,R0,V0);
              _TBD.CalcuSfAccCentral3(ord/2+1);
            }
        }while(CmpFT(ord/2)>epslong);

      _TBD.CalcuCentralSf3(ord,R0,V0);
      _TBD.CalcuSfAccCentral3(ord/2+1);

      _TBD.CalcuIntXv(ord/2+1,ord,ord/2+1,R,V);   A=CalcuAcc(R);
      _TBD.WriteF(ord/2+1,A);
      _TBD.CalcuIntXv(ord/2+1,ord,ord/2,R,V);   A=CalcuAcc(R);
      _TBD.WriteF(ord/2+1,A);

      _TBD.CalcuIntXv(-ord/2-1,ord,-ord/2-1,R,V);   A=CalcuAcc(R);
      _TBD.WriteF(-ord/2-1,A);
      _TBD.CalcuIntXv(-ord/2-1,ord,-ord/2,R,V);   A=CalcuAcc(R);
      _TBD.WriteF(-ord/2-1,A);
    }

  _TBD.CalcuCentralSf3(order,R0,V0);
  _TBD.CalcuSfAccCentral3(7); 
}

void sbs_ord(int ord)
{ TVector R,V,A;

  _TBD.sbs(1,7);
  _TBD.CalcuIntXv(7,ord,ord/2+1,R,V);
  A=CalcuAcc(R);
  _TBD.WriteF(7,A);

  _TBD.MoveForwardOnce();
  _TC+=_TBD.h;      

  _TBD.CalcuIntXv(6,ord,ord/2,R,V);
  A=CalcuAcc(R);
  _TBD.WriteF(6,A);

}


double FindMidValue(double v[])
{ for(int i=0;i<5;i++)
    { int c=0;
      for(int j=0;j<5;j++)
         if(v[j]>v[i]) c++;
      if(c==2) return v[i];
    }
  return v[0];
}



int main()
{ _mu=SunGm;
  _InitialElem[0]=1;  _InitialElem[1]=0.01; for(int i=2;i<6;i++) _InitialElem[i]=1;
  //-*+*+*+*+*+**+*+*-*+*+*+*+*+**+*+*-*+*+*+*+*+**+*+*-*+*+*+*+*+**+*+*-*+*+*+

  ofstream out("ddd.txt");  out<<setprecision(18);
  TVector R0,V0,R,V,Ra,Va;
  double h,max,elem[6],b[5],dr,mid,p,pb,Period;
  int cnt;

  Period=6.283185307179586476925286766559L/sqrtl(_mu/(_InitialElem[0]*_InitialElem[0]*_InitialElem[0]));

  h=1;
  Initia(12,h);

  max=-1;  cnt=0; pb=0;
  while(max<0.01)
    { _TBD.CalcuIntXv(0,12,0,R,V);
      ElemToCor(_mu,_TC,_InitialElem,Ra,Va);

      dr=R.DistanceTo(Ra);
      b[cnt++]=dr;  cnt%=5;

      if(_TC>h*6)
        { mid=FindMidValue(b);
          if(max<mid) max=mid;

          p=log10(_TC/Period);
          if(p>pb)
            { pb+=0.2;

              CorToElem(_mu,-_TC,R,V,elem);  for(int i=0;i<6;i++) elem[i]-=_InitialElem[i];
              out<<p<<" "<<_TC<<" "<<max;
              for(int i=0;i<6;i++) out<<" "<<elem[i];
              out<<endl;
            }
        }
      sbs_ord(12);
    }
  return 0;
}
//---------------------------------------------------------------------------




