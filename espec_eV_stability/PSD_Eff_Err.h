
class PSD_Eff_Err{
  double psd_er_d, psd_er_u;
 public:
  
  double psd_eff,n_total_tag;
  double PSD_ER_D()
  {
    if(psd_eff==0.0) 
      { 
	psd_er_d = 0.0; 
      }
    
    if(psd_eff==1.0)
      { 
	psd_er_d = 1.0-(0.32**((1.0/(n_total_tag+1.0)))); 
      }
    
    if(psd_eff>0.0&&psd_eff<1.0)
      { 
	psd_er_d = (1.0/(sqrt(n_total_tag)))*sqrt((psd_eff*(1.0-psd_eff)));
	
      }
    return psd_er_d; 
  }
  
double PSD_ER_U()
  {
    if(psd_eff==0.0) 
      { 
	psd_er_u = 1.0-(0.32**((1.0/(n_total_tag+1.0)))); 
      }
    
    if(psd_eff==1.0)
      { 
	psd_er_u = 0.0;
      }
    
    if(psd_eff>0.0&&psd_eff<1.0)
      { 
	psd_er_u = (1.0/(sqrt(n_total_tag)))*sqrt((psd_eff*(1.0-psd_eff)));
	
      }
    return psd_er_u; 
  }

};
