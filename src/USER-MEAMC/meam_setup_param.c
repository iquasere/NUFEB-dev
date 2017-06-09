#include "meam.h"

//
//     do a sanity check on index parameters
      void meam_checkindex(int num, int lim, int nidx, int *idx /*idx(3)*/, int *ierr)
      {
        //: idx[0..2]
        *ierr = 0;
        if (nidx < num) {
           *ierr = 2;
           return;
        }

        for (int i=0; i<num; i++) {
          if ((idx[i] < 1) || (idx[i] > lim)) {
              *ierr = 3;
              return;
          }
        }
      }

//
//     Declaration in pair_meam.h:
//
//     void meam_setup_param(int *, double *, int *, int *, int *);
//
//     in pair_meam.cpp
//
//     meam_setup_param(&which,&value,&nindex,index,&errorflag);
//
//
//
//     The "which" argument corresponds to the index of the "keyword" array
//     in pair_meam.cpp:
//
//     0 = Ec_meam
//     1 = alpha_meam
//     2 = rho0_meam
//     3 = delta_meam
//     4 = lattce_meam
//     5 = attrac_meam
//     6 = repuls_meam
//     7 = nn2_meam
//     8 = Cmin_meam
//     9 = Cmax_meam
//     10 = rc_meam
//     11 = delr_meam
//     12 = augt1
//     13 = gsmooth_factor
//     14 = re_meam
//     15 = ialloy
//     16 = mixture_ref_t
//     17 = erose_form
//     18 = zbl_meam
//     19 = emb_lin_neg
//     20 = bkgd_dyn

      void meam_setup_param_(int *which_p, double *value_p, int *nindex_p, int *index /*index(3)*/, int *errorflag)
      {
        //: index[0..2]
        int i1, i2;
        *errorflag = 0;
        int which = *which_p;
        double value = *value_p;
        int nindex = *nindex_p;

      switch(which) {
  //     0 = Ec_meam
        case 0:
          meam_checkindex(2,maxelt,nindex,index,errorflag);
          if (*errorflag!=0) return;
          meam_data.Ec_meam[index[0]][index[1]] = value;
          break;

  //     1 = alpha_meam
        case 1:
          meam_checkindex(2,maxelt,nindex,index,errorflag);
          if (*errorflag!=0) return;
          meam_data.alpha_meam[index[0]][index[1]] = value;
          break;

  //     2 = rho0_meam
        case 2:
          meam_checkindex(1,maxelt,nindex,index,errorflag);
          if (*errorflag!=0) return;
          meam_data.rho0_meam[index[0]] = value;
          break;

  //     3 = delta_meam
        case 3:
          meam_checkindex(2,maxelt,nindex,index,errorflag);
          if (*errorflag!=0) return;
          meam_data.delta_meam[index[0]][index[1]] = value;
          break;

  //     4 = lattce_meam
        case 4:
          meam_checkindex(2,maxelt,nindex,index,errorflag);
          if (*errorflag!=0) return;
          int val = (int)value;

          if (val==0)
            meam_data.lattce_meam[index[0]][index[1]] = FCC;
          else if (val==1)
            meam_data.lattce_meam[index[0]][index[1]] = BCC;
          else if (val==2)
            meam_data.lattce_meam[index[0]][index[1]] = HCP;
          else if (val==3)
            meam_data.lattce_meam[index[0]][index[1]] = DIM;
          else if (val==4)
            meam_data.lattce_meam[index[0]][index[1]] = DIA;
          else if (val==5)
            meam_data.lattce_meam[index[0]][index[1]] = B1;
          else if (val==6)
            meam_data.lattce_meam[index[0]][index[1]] = C11;
          else if (val==7)
            meam_data.lattce_meam[index[0]][index[1]] = L12;
          else if (val==8)
            meam_data.lattce_meam[index[0]][index[1]] = B2;
          break;

  //     5 = attrac_meam
        case 5:
          meam_checkindex(2,maxelt,nindex,index,errorflag);
          if (*errorflag!=0) return;
          meam_data.attrac_meam[index[0]][index[1]] = value;
          break;

  //     6 = repuls_meam
        case 6:
          meam_checkindex(2,maxelt,nindex,index,errorflag);
          if (*errorflag!=0) return;
          meam_data.repuls_meam[index[0]][index[1]] = value;
          break;

  //     7 = nn2_meam
        case 7:
          meam_checkindex(2,maxelt,nindex,index,errorflag);
          if (*errorflag!=0) return;
          i1 = min(index[0],index[1]);
          i2 = max(index[0],index[1]);
          meam_data.nn2_meam[i1][i2] = (int)value;
          break;

  //     8 = Cmin_meam
        case 8:
          meam_checkindex(3,maxelt,nindex,index,errorflag);
          if (*errorflag!=0) return;
          meam_data.Cmin_meam[index[0]][index[1]][index[2]] = value;
          break;

  //     9 = Cmax_meam
        case 9:
          meam_checkindex(3,maxelt,nindex,index,errorflag);
          if (*errorflag!=0) return;
          meam_data.Cmax_meam[index[0]][index[1]][index[2]] = value;
          break;

  //     10 = rc_meam
        case 10:
          meam_data.rc_meam = value;
          break;

  //     11 = delr_meam
        case 11:
          meam_data.delr_meam = value;
          break;

  //     12 = augt1
        case 12:
          meam_data.augt1 = (int)value;
          break;

  //     13 = gsmooth
        case 13:
          meam_data.gsmooth_factor = value;
          break;

  //     14 = re_meam
        case 14:
          meam_checkindex(2,maxelt,nindex,index,errorflag);
          if (*errorflag!=0) return;
          meam_data.re_meam[index[0]][index[1]] = value;
          break;

  //     15 = ialloy
        case 15:
          meam_data.ialloy = (int)value;
          break;

  //     16 = mixture_ref_t
        case 16:
          meam_data.mix_ref_t = (int)value;
          break;

  //     17 = erose_form
        case 17:
          meam_data.erose_form = (int)value;
          break;

  //     18 = zbl_meam
        case 18:
          meam_checkindex(2,maxelt,nindex,index,errorflag);
          if (*errorflag!=0) return;
          i1 = min(index[0],index[1]);
          i2 = max(index[0],index[1]);
          meam_data.zbl_meam[i1][i2] = (int)value;
          break;

  //     19 = emb_lin_neg
        case 19:
          meam_data.emb_lin_neg = (int)value;
          break;

  //     20 = bkgd_dyn
        case 20:
          meam_data.bkgd_dyn = (int)value;
          break;

        default:
          *errorflag = 1;
      }
    }
