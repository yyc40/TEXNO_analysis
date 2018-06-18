double ERR_CAL(double S1,double er_S1,int N1,double er_N1,double S2,double er_S2,int N2,double er_N2) 
{
  double er_cal=0;
  double er_cal1=0;
  double er_cal2=0;
  if (N1!=0) {
    er_cal1=(S1*(double)N1)*sqrt(((er_S1/S1)*(er_S1/S1))+((er_N1/(double)N1)*(er_N1/(double)N1)));
  }
  else {
    er_cal1=S1;
  }
  if (N2!=0) {
    er_cal2=(S2*(double)N2)*sqrt(((er_S2/S2)*(er_S2/S2))+((er_N2/(double)N2)*(er_N2/(double)N2)));
  }
  else {
    er_cal2=S2;
  }
  
  er_cal=sqrt(er_cal1*er_cal1+er_cal2*er_cal2);
  return er_cal;
}


double Total_ERR_U(double S,double er_S,double Eff,double er_Eff) 
{
  double er_cal_u;
  if(S!=0&&Eff!=0) 
    {
      er_cal_u=(S/Eff)*sqrt(((er_S/S)*(er_S/S))+((er_Eff/Eff)*(er_Eff/Eff)));
    }
  if(S==0.0)
    {
      er_cal_u=Eff;
    }

  return er_cal_u;
}

double Total_ERR_D(double S,double er_S,double Eff,double er_Eff) 
{
  double er_cal_d;
  if(S!=0&&Eff!=0) 
    {
      er_cal_d=(S/Eff)*sqrt(((er_S/S)*(er_S/S))+((er_Eff/Eff)*(er_Eff/Eff)));
    }
  if(S==0.0)
    {
      er_cal_d=0.0;
    }

  return er_cal_d;
}
