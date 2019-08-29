
#include <RcppEigen.h>
#include <Eigen/Core>
#include <iostream.h>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatD;
// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXd invvech_eigen(const Eigen::VectorXd & x){
  int k=0;
  int leng = x.size();
  int tam=(sqrt(8*leng+1)-1)/2;
  Eigen::MatrixXd mat(tam, tam);
  for(int i=0; i<tam; ++i){
    for(int j=i; j<tam; ++j){
      mat(i,j)=x[k];
      k=k+1;
    }
  }
  return mat.selfadjointView<Eigen::Upper>();
}

//https://stackoverflow.com/questions/28950857/how-to-construct-block-diagonal-matrix
//Construção de matriz bloco diagonal somente com uma matriz
// [[Rcpp::export]]
Eigen::MatrixXd blkdiag(const Eigen::MatrixXd & a, int count){
  Eigen::MatrixXd bdm = Eigen::MatrixXd::Zero(a.rows() * count, a.cols() * count);
  for (int i = 0; i < count; ++i){
    bdm.block(i * a.rows(), i * a.cols(), a.rows(), a.cols()) = a;
  }
  return bdm;
}

// [[Rcpp::export]]
Rcpp::List initial_beta(
    int nvar,
    const Eigen::MatrixXd & x,
    const Eigen::VectorXd & y,
    int ncluster,
    const Eigen::MatrixXd & panelsetup2,
    const Eigen::VectorXd & wi_j_star,
    const Eigen::VectorXd & wj_star){

  Eigen::MatrixXd somat1 = Eigen::MatrixXd::Zero(nvar,nvar);
  Eigen::VectorXd somat3 = Eigen::VectorXd::Zero(nvar);

  for(int i = 0; i < ncluster; ++i){
    int b = panelsetup2(i,0);
    int a = panelsetup2(i,1);
    Eigen::VectorXd yj = y.segment(a,b);
    Eigen::MatrixXd xj = x.block(a,0,b,nvar);

    Eigen::VectorXd v = wi_j_star.segment(a,b);
    Eigen::MatrixXd diag_ = v.asDiagonal();
    double wj = wj_star(i);

    Eigen::MatrixXd wj_t_xj_diag = wj*xj.transpose()*diag_;
    Eigen::MatrixXd t1 = wj_t_xj_diag*xj;
    somat1 += t1;

    Eigen::VectorXd t3 = wj_t_xj_diag*yj;
    somat3 += t3;
  }

  return Rcpp::List::create(Rcpp::Named("somat1")=somat1,
                            Rcpp::Named("somat3")=somat3);
}

// [[Rcpp::export]]
Rcpp::List iterative_multi_beta(
              int nvar,
              int ncluster,
              const Eigen::MatrixXd & panelsetup2,
              const Eigen::MatrixXd & x,
              const Eigen::VectorXd & y,
              const Eigen::VectorXd & z,
              const Eigen::VectorXd & wi_j_star,
              const Eigen::VectorXd & wj_star,
              const Eigen::MatrixXd & teta,
              int tt,
              const Eigen::MatrixXd & teta_genlin){

  Eigen::MatrixXd s_matp = Eigen::MatrixXd::Zero(nvar,nvar);
  Eigen::VectorXd s_matq = Eigen::VectorXd::Zero(nvar);
  std::vector<Eigen::MatrixXd> invvjs;
  Eigen::MatrixXd sigma_inv = teta_genlin.inverse();
  Eigen::MatrixXd s_sigma2e = teta.block(0,0,1,1).inverse();

  for(int i = 0; i < ncluster; ++i){
    int b = panelsetup2(i,0);
    int a = panelsetup2(i,1);
    Eigen::VectorXd yj = y.segment(a,b);
    Eigen::MatrixXd xj = x.block(a,0,b,nvar);
    Eigen::VectorXd zj = z.segment(a,b);

    Eigen::VectorXd v = wi_j_star.segment(a,b);
    Eigen::MatrixXd diag_ = v.asDiagonal();

    double wj = wj_star(i);

    int np_div_tt = b/tt;
    Eigen::MatrixXd diag_solve_sigma = diag_*blkdiag(sigma_inv, np_div_tt);
    Eigen::VectorXd diag_solve_sigma_zj = diag_solve_sigma*zj;
    Eigen::MatrixXd aj = (s_sigma2e+zj.transpose()*diag_solve_sigma_zj).inverse();
    Eigen::MatrixXd invvj = diag_solve_sigma - diag_solve_sigma_zj*aj*diag_solve_sigma_zj.transpose();

    invvjs.push_back (invvj);

    Eigen::MatrixXd wj_t_xj_invvj = wj*xj.transpose()*invvj;
    Eigen::MatrixXd t1 = wj_t_xj_invvj*xj;
    s_matp += t1;

    Eigen::VectorXd t2 = wj_t_xj_invvj*yj;
    s_matq += t2;
  }

  return Rcpp::List::create(Rcpp::Named("s_matp")=s_matp,
                            Rcpp::Named("s_matq")=s_matq,
                            Rcpp::Named("invvjs")=invvjs);
}

// [[Rcpp::export]]
Rcpp::List iterative_multi_teta(
    int ncluster,
    const Eigen::VectorXd & beta,
    const Eigen::MatrixXd & panelsetup2,
    const Eigen::MatrixXd & x,
    const Eigen::VectorXd & y,
    const Eigen::VectorXd & wi_j_star,
    const Eigen::VectorXd & inv_wi_j_star,
    const Eigen::VectorXd & wj_star,
    int tt,
    const std::vector<Eigen::MatrixXd> & DELTA,
    const std::vector<Eigen::MatrixXd> & invvjs,
    const Eigen::VectorXd & h_matrix,
    int s,
    int ss,
    int nvar){

  Eigen::MatrixXd R(ncluster,ss);
  Eigen::MatrixXd S(ncluster,s);

  for(int i = 0; i < ncluster; ++i){
    Eigen::MatrixXd Rklj(s,s);
    Eigen::VectorXd Skj(s);

    int b = panelsetup2(i,0);
    int a = panelsetup2(i,1);
    Eigen::VectorXd yj = y.segment(a,b); //Outro formato de matriz 'info_j' era melhor
    Eigen::MatrixXd xj = x.block(a,0,b,nvar);

    Eigen::VectorXd v = wi_j_star.segment(a,b);
    Eigen::MatrixXd diag_ = v.asDiagonal();

    double wj = wj_star(i);
    Eigen::MatrixXd invvj = invvjs[i];

    Eigen::MatrixXd zj_zj = Eigen::MatrixXd::Ones(b,b);
    int np_div_tt = b/tt;
    Eigen::VectorXd inv_v = inv_wi_j_star.segment(a,b);
    Eigen::MatrixXd s_diag = inv_v.asDiagonal();

    Eigen::VectorXd eij = yj - xj*beta;

   for(int k = 0; k < s; ++k){
     Eigen::MatrixXd Rklj1 = invvj*(h_matrix[k]*zj_zj + blkdiag(DELTA[k], np_div_tt)*s_diag)*invvj;
     Skj(k)= wj*(Rklj1*(eij*eij.transpose())).trace();
     for(int l = 0; l < s; ++l){
       Rklj(k,l)= wj*(Rklj1*(h_matrix[l]*zj_zj + blkdiag(DELTA[l], np_div_tt)*s_diag)).trace();
     }
   }

    Eigen::VectorXd t11(Eigen::Map<Eigen::VectorXd>(Rklj.data(), ss));

    R.row(i) = t11;
    S.row(i) = Skj;
  }

  return Rcpp::List::create(Rcpp::Named("R")=R,
                            Rcpp::Named("S")=S);
}

// [[Rcpp::export]]
Rcpp::List multi_variances_residuals(
    int nvar,
    int s,
    const Eigen::VectorXd beta,
    const Eigen::VectorXd teta,
    const Eigen::VectorXd & y,
    const Eigen::MatrixXd & x,
    const Eigen::VectorXd & z,
    const Eigen::VectorXd & wi_j_star,
    int ncluster,
    const Eigen::MatrixXd & panelsetup2,
    const Eigen::MatrixXd & R,
    const Eigen::MatrixXd & S,
    const Eigen::VectorXd & wj_star,
    int tt,
    const Eigen::MatrixXd & teta_1,
    const Eigen::MatrixXd & teta_matrix,
    int nsubjc_t,
    double teta1){

  Eigen::MatrixXd s_mat_c = Eigen::MatrixXd::Zero(nvar,nvar);
  Eigen::MatrixXd s_mat_d = Eigen::MatrixXd::Zero(s,s);

  Eigen::VectorXd u(ncluster);
  Eigen::VectorXd var_u(ncluster);
  Eigen::VectorXd yhat(nsubjc_t);
  Eigen::VectorXd v(nsubjc_t);

  for(int i = 0; i < ncluster; ++i){
    int b = panelsetup2(i,0);
    int a = panelsetup2(i,1);
    Eigen::VectorXd yj = y.segment(a,b); //Outro formato de matriz 'info_j' era melhor
    Eigen::MatrixXd xj = x.block(a,0,b,nvar);
    Eigen::VectorXd zj = z.segment(a,b);

    double wj = wj_star(i);
    double wj_quad = std::pow(wj,2);

    Eigen::VectorXd v2 = wi_j_star.segment(a,b);
    Eigen::MatrixXd diag_ = v2.asDiagonal();

    Eigen::VectorXd xj_beta = xj*beta;
    Eigen::VectorXd eij = yj - xj_beta;

    int np_div_tt = b/tt;
    Eigen::MatrixXd sigma = blkdiag(teta_matrix, np_div_tt);

    Eigen::MatrixXd diag_sigma_inverse = diag_*sigma.inverse();
    Eigen::MatrixXd aj = (teta_1 + zj.transpose()*diag_sigma_inverse*zj).inverse();
    Eigen::MatrixXd invvj = diag_sigma_inverse - diag_sigma_inverse*zj*aj*zj.transpose()*diag_sigma_inverse;

    Eigen::VectorXd cj = xj.transpose()*invvj*eij;
    Eigen::RowVectorXd t_cj = cj;
    Eigen::MatrixXd wj_quad_cj = wj_quad*(cj*t_cj);
    s_mat_c += wj_quad_cj;

    Eigen::VectorXd R_row = R.row(i);
    Eigen::MatrixXd Rklj(Eigen::Map<Eigen::MatrixXd>(R_row.data(), s, s));
    Eigen::VectorXd Skj = S.row(i);

    Eigen::VectorXd Rklj_teta = Rklj*teta;
    Eigen::VectorXd SRt = Skj - Rklj_teta;

    Eigen::MatrixXd Dkj = SRt*SRt.transpose();
    s_mat_d += Dkj;

    //Residuals

    Eigen::RowVectorXd Rhj = Eigen::VectorXd::Constant(b, teta1);

    Eigen::MatrixXd Vj = zj*Rhj + sigma;
    Eigen::VectorXd aux = Rhj*Vj.inverse();

    double aux_eij = aux.dot(eij);
    u(i) = aux_eij;

    double aux_rhj = aux.dot(Rhj);
    double diag_var_u = teta1 - aux_rhj;
    var_u(i) = diag_var_u;

    Eigen::VectorXd zj_aux_eij = zj*aux_eij; //aux_eij é um double
    Eigen::VectorXd yhat1 = xj_beta + zj_aux_eij;
    yhat.segment(a,b) = yhat1;

    Eigen::VectorXd vj = eij - zj_aux_eij;
    v.segment(a,b) = vj;
  }
  return Rcpp::List::create(Rcpp::Named("s_mat_c")=s_mat_c,
                            Rcpp::Named("s_mat_d")=s_mat_d,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("var_u")=var_u,
                            Rcpp::Named("v")=v,
                            Rcpp::Named("yhat")=yhat);
}

// [[Rcpp::export]]
Rcpp::List tjs_uni_beta(int p,
                        int m,
                        const Eigen::MatrixXd & panelsetup,
                        const Eigen::MatrixXd & x,
                        const Eigen::VectorXd & y,
                        const Eigen::MatrixXd & z,
                        const Eigen::VectorXd & wi_j_star,
                        int q,
                        int s,
                        const Eigen::MatrixXd & H,
                        const Eigen::VectorXd & wj_star){
  std::vector<Eigen::MatrixXd> T1;
  std::vector<Eigen::MatrixXd> T2;
  std::vector<Eigen::VectorXd> T3;
  std::vector<Eigen::VectorXd> T4;
  std::vector<Eigen::MatrixXd> T5;
  std::vector<Eigen::MatrixXd> H_K;
  std::vector<Eigen::MatrixXd> T5_HK;
  std::vector<double> TR_T5_HK;
  Eigen::MatrixXd somat1 = Eigen::MatrixXd::Zero(p,p);
  Eigen::VectorXd somat3 = Eigen::VectorXd::Zero(p);

  for(int k = 0; k < s; ++k){
    Eigen::VectorXd h_k = H.row(k);
    Eigen::MatrixXd H_k(Eigen::Map<Eigen::MatrixXd>(h_k.data(), q, q));

    H_K.push_back (H_k);
  }

  for(int j = 0; j < m; ++j){
    int b = panelsetup(j,0);
    int a = panelsetup(j,1);
    Eigen::VectorXd yj = y.segment(a,b);
    Eigen::MatrixXd xj = x.block(a,0,b,p);
    Eigen::MatrixXd zj = z.block(a,0,b,q);

    Eigen::VectorXd v = wi_j_star.segment(a,b);
    Eigen::MatrixXd diag = v.asDiagonal();

    Eigen::MatrixXd t_xj_diag = xj.transpose()*diag;
    Eigen::MatrixXd t_zj_diag = zj.transpose()*diag;

    Eigen::MatrixXd t1j = t_xj_diag*xj;
    Eigen::MatrixXd t2j = t_xj_diag*zj;
    Eigen::VectorXd t3j = t_xj_diag*yj;
    Eigen::VectorXd t4j = t_zj_diag*yj;
    Eigen::MatrixXd t5j = t_zj_diag*zj;

    for(int l = 0; l < s; ++l){
      Eigen::MatrixXd t5j_hk = t5j*H_K[l];
      T5_HK.push_back (t5j_hk);

      double tr_t5j_hk = t5j_hk.trace();
      TR_T5_HK.push_back (tr_t5j_hk);
    }
    //initial beta
    double wj = wj_star(j);

    Eigen::MatrixXd t11 = wj*t1j;
    somat1 += t11;

    Eigen::VectorXd t33 = wj*t3j;
    somat3 += t33;

    T1.push_back (t1j);
    T2.push_back (t2j);
    T3.push_back (t3j);
    T4.push_back (t4j);
    T5.push_back (t5j);
  }

  return Rcpp::List::create(Rcpp::Named("T1")=T1,
                            Rcpp::Named("T2")=T2,
                            Rcpp::Named("T3")=T3,
                            Rcpp::Named("T4")=T4,
                            Rcpp::Named("T5")=T5,
                            Rcpp::Named("H_K")=H_K,
                            Rcpp::Named("T5_HK")=T5_HK,
                            Rcpp::Named("TR_T5_HK")=TR_T5_HK,
                            Rcpp::Named("somat1")=somat1,
                            Rcpp::Named("somat3")=somat3);
}

// [[Rcpp::export]]
Rcpp::List iterative_uni_beta(int nvar,
                               int ncluster,
                               const Eigen::MatrixXd & tsit,
                               const std::vector<Eigen::MatrixXd> & T1,
                               const std::vector<Eigen::MatrixXd> & T2,
                               const std::vector<Eigen::VectorXd> & T3,
                               const std::vector<Eigen::VectorXd> & T4,
                               const std::vector<Eigen::MatrixXd> & T5,
                               const Eigen::VectorXd & wj_star){
  Eigen::MatrixXd s_matp = Eigen::MatrixXd::Zero(nvar,nvar);
  Eigen::VectorXd s_matq = Eigen::VectorXd::Zero(nvar);
  std::vector<Eigen::MatrixXd> AJS;
  for(int j = 0; j < ncluster; ++j){
    double wj = wj_star(j);

    Eigen::MatrixXd Aj = (T5[j] + tsit).inverse();
    AJS.push_back (Aj);

    Eigen::MatrixXd t0 = T2[j]*Aj;
    Eigen::MatrixXd t1 = wj*(T1[j] - t0*T2[j].transpose());
    s_matp += t1;
    Eigen::VectorXd t2 = wj*(T3[j] - t0*T4[j]);
    s_matq += t2;
  }

  return Rcpp::List::create(Rcpp::Named("AJS")=AJS,
                            Rcpp::Named("s_matp")=s_matp,
                            Rcpp::Named("s_matq")=s_matq);
}

// [[Rcpp::export]]
Rcpp::List iterative_uni_theta(int s,
                               int nvar,
                               int ncluster,
                               const Eigen::MatrixXd & sit,
                               const Eigen::MatrixXd & tsit,
                               const Eigen::MatrixXd & panelsetup,
                               const Eigen::VectorXd & y,
                               const Eigen::MatrixXd & x,
                               const Eigen::MatrixXd & z,
                               int q,
                               const Eigen::VectorXd & beta,
                               const Eigen::VectorXd & wj_star,
                               double teta_s,
                               const std::vector<Eigen::MatrixXd> & AJS,
                               const std::vector<Eigen::MatrixXd> & T5,
                               const Eigen::VectorXd & wi_j_star,
                               const std::vector<Eigen::MatrixXd> & H_K,
                               const std::vector<double> & TR_T5_HK,
                               const std::vector<Eigen::MatrixXd> & T5_HK){

  Eigen::MatrixXd r_mat = Eigen::MatrixXd::Zero(s,s);
  Eigen::VectorXd s_mat = Eigen::VectorXd::Zero(s);

  std::vector<int> delta(s, 0);
  delta.back() = 1;

  for(int j = 0; j < ncluster; ++j){
    int b = panelsetup(j,0);
    int a = panelsetup(j,1);
    Eigen::VectorXd yj = y.segment(a,b); //Outro formato de matriz 'info_j' era melhor
    Eigen::MatrixXd xj = x.block(a,0,b,nvar);
    Eigen::MatrixXd zj = z.block(a,0,b,q);

    double wj = wj_star(j);

    Eigen::VectorXd eij = yj - xj*beta;
    Eigen::RowVectorXd t_eij = eij;

    Eigen::MatrixXd Rklj(s, s);
    Eigen::VectorXd Skj(s);

    Eigen::MatrixXd aj = AJS[j];
    Eigen::MatrixXd t5j = T5[j];
    Eigen::MatrixXd teta_s_aj_sit = teta_s*aj*sit;
    Eigen::MatrixXd t5j_aj = t5j*aj;

    Eigen::VectorXd v = wi_j_star.segment(a,b);
    Eigen::MatrixXd DIAG = v.asDiagonal();

    Eigen::MatrixXd t_zj_diag = zj.transpose()*DIAG;
    double tr_e_diag_e = (t_eij*DIAG*eij).trace();
    Eigen::VectorXd zj_diag_e = t_zj_diag*eij;
    Eigen::RowVectorXd t_zj_diag_e = zj_diag_e;

    for(int k = 0; k < s; ++k){
      int cont = j*s;
      int delta_k = delta[k];

      Eigen::MatrixXd B_k = teta_s_aj_sit*H_K[k] - delta_k*aj;
      Eigen::MatrixXd C_k = -delta_k*aj + B_k - B_k*t5j_aj;

      Skj(k) = wj*(delta_k*tr_e_diag_e +(t_zj_diag_e*C_k*zj_diag_e).trace());

      Eigen::MatrixXd t5j_Ck = t5j*C_k;
      double tr_t5j_Ck = t5j_Ck.trace();

      for(int l = 0; l < s; ++l){
        int delta_l = delta[l];
        Rklj(k,l) = wj*(delta_k*delta_l*b + delta_l*tr_t5j_Ck + delta_k*TR_T5_HK[cont] + (t5j_Ck*T5_HK[cont]).trace());
        cont += 1;
      }
    }
    r_mat += Rklj;
    s_mat += Skj;
  }

  return Rcpp::List::create(Rcpp::Named("r_mat")=r_mat,
                            Rcpp::Named("s_mat")=s_mat);
}

// [[Rcpp::export]]
Rcpp::List scaled_weight(
    int n,
    int m,
    const Eigen::MatrixXd & panelsetup,
    const Eigen::VectorXd & wi_j,
    const Eigen::VectorXd & wj){

  Eigen::VectorXd wi_j_star(n);
  Eigen::VectorXd cluster_wgt(m);

  for(int j = 0; j < m; ++j){
    int b = panelsetup(j,0);
    int a = panelsetup(j,1);

    wi_j_star.segment(a,b) = wi_j.segment(a,b)/(wi_j.segment(a,b)).mean();
    cluster_wgt(j) = wj(a);
  }

  Eigen::VectorXd wj_star = cluster_wgt/cluster_wgt.mean();

  return Rcpp::List::create(Rcpp::Named("wi_j_star")=wi_j_star,
                            Rcpp::Named("wj_star")=wj_star);
}

// [[Rcpp::export]]
Rcpp::List initial_theta(
    const Eigen::VectorXd & beta0,
    const Eigen::MatrixXd & x,
    const Eigen::VectorXd & y,
    int ncluster,
    const Eigen::MatrixXd & panelsetup,
    const Eigen::VectorXd & wi_j_star,
    const Eigen::VectorXd & wj_star,
    int nvar){

  double wj_t6 = 0;
  double wj_aux = 0;

  for(int j = 0; j < ncluster; ++j){
    int b = panelsetup(j,0);
    int a = panelsetup(j,1);
    Eigen::VectorXd yj = y.segment(a,b); //Outro formato de matriz 'info_j' era melhor
    Eigen::MatrixXd xj = x.block(a,0,b,nvar);

    Eigen::RowVectorXd wi_j_starj = wi_j_star.segment(a,b);
    double wj = wj_star(j);

    Eigen::VectorXd resid = yj - xj*beta0;
    Eigen::MatrixXd ebarj = (wi_j_starj*resid)/b;
    double k = ebarj.value();
    Eigen::VectorXd resbar(b);
    for(int nj = 0; nj < b; ++nj){
      resbar(nj) = std::pow(resid(nj) - k, 2);
    }
    wj_t6 += wj*((wi_j_starj*resbar).value());
    wj_aux += wj*(b - 1);
  }

  return Rcpp::List::create(Rcpp::Named("wj_t6")=wj_t6,
                            Rcpp::Named("wj_aux")=wj_aux);
}

// [[Rcpp::export]]
Rcpp::List uni_variances_residuals(
    int nvar,
    int s,
    const Eigen::VectorXd beta,
    const Eigen::VectorXd teta,
    const Eigen::VectorXd & y,
    const Eigen::MatrixXd & x,
    int ncluster,
    const Eigen::MatrixXd & panelsetup,
    const std::vector<Eigen::MatrixXd> & T2,
    const std::vector<Eigen::MatrixXd> & T5,
    const Eigen::VectorXd & wj_star,
    const Eigen::MatrixXd & sit,
    int nsubjc_t,
    int q,
    const Eigen::MatrixXd & z,
    const Eigen::MatrixXd & Sigmau,
    double teta_s,
    const Eigen::VectorXd & wi_j_star,
    const std::vector<Eigen::MatrixXd> & H_K,
    const std::vector<double> & TR_T5_HK,
    const std::vector<Eigen::MatrixXd> & T5_HK){

  Eigen::MatrixXd u(ncluster, q);
  Eigen::MatrixXd var_u(ncluster, q);
  Eigen::VectorXd yhat(nsubjc_t);
  Eigen::VectorXd v(nsubjc_t);
  Eigen::VectorXd var_v(nsubjc_t);

  Eigen::MatrixXd s_matc = Eigen::MatrixXd::Zero(nvar, nvar);
  Eigen::MatrixXd s_matd = Eigen::MatrixXd::Zero(s, s);

  Eigen::MatrixXd tsit = teta_s*sit;
  std::vector<int> delta(s, 0);//-
  delta.back() = 1;//-

  for(int j = 0; j < ncluster; ++j){
    int b = panelsetup(j,0);
    int a = panelsetup(j,1);
    Eigen::VectorXd yj = y.segment(a,b); //Outro formato de matriz 'info_j' era melhor
    Eigen::MatrixXd xj = x.block(a,0,b,nvar);
    Eigen::MatrixXd zj = z.block(a,0,b,q);

    //Variances
    double wj = wj_star(j);
    double wj_quad = std::pow(wj,2);

    Eigen::VectorXd xj_beta = xj*beta;
    Eigen::VectorXd eij = yj - xj_beta;
    Eigen::RowVectorXd t_eij = eij;//-

    Eigen::VectorXd v2 = wi_j_star.segment(a,b);//-
    Eigen::MatrixXd DIAG = v2.asDiagonal();//-

    Eigen::VectorXd t7j = xj.transpose()*DIAG*eij;
    Eigen::MatrixXd t_zj_diag = zj.transpose()*DIAG;
    Eigen::VectorXd t8j = t_zj_diag*eij;

    Eigen::MatrixXd t5j = T5[j];//-
    Eigen::MatrixXd aj = (t5j + tsit).inverse();//-

    Eigen::MatrixXd teta_s_aj_sit = aj*tsit;//-
    Eigen::MatrixXd t5j_aj = t5j*aj;//-
    Eigen::VectorXd cj = t7j - T2[j]*aj*t8j;

    Eigen::MatrixXd wj_cj = wj_quad*(cj*cj.transpose());
    s_matc += wj_cj;

    double tr_e_diag_e = (t_eij*DIAG*eij).trace();//-
    Eigen::VectorXd zj_diag_e = t_zj_diag*eij;//-
    Eigen::RowVectorXd t_zj_diag_e = zj_diag_e;//-

    Eigen::MatrixXd Rklj(s, s);
    Eigen::VectorXd Skj(s);

    for(int k = 0; k < s; ++k){
      int cont = j*s;
      int delta_k = delta[k];

      Eigen::MatrixXd B_k = teta_s_aj_sit*H_K[k] - delta_k*aj;
      Eigen::MatrixXd C_k = -delta_k*aj + B_k - B_k*t5j_aj;

      Skj(k) = wj*(delta_k*tr_e_diag_e +(t_zj_diag_e*C_k*zj_diag_e).trace());

      Eigen::MatrixXd t5j_Ck = t5j*C_k;
      double tr_t5j_Ck = t5j_Ck.trace();

      for(int l = 0; l < s; ++l){
        int delta_l = delta[l];
        Rklj(k,l) = wj*(delta_k*delta_l*b + delta_l*tr_t5j_Ck + delta_k*TR_T5_HK[cont] + (t5j_Ck*T5_HK[cont]).trace());
        cont += 1;
      }
    }

    Eigen::VectorXd Rklj_teta = Rklj*teta;
    Eigen::VectorXd SRt = Skj - Rklj_teta;

    Eigen::MatrixXd Dkj = wj_quad*(SRt*SRt.transpose());
    s_matd += Dkj;

    //Residuals
    Eigen::MatrixXd t_zj = zj.transpose();

    Eigen::MatrixXd Rhj = Sigmau*t_zj;
    Eigen::MatrixXd t_Rhj = Rhj.transpose();
    Eigen::MatrixXd Vj = t_Rhj*t_zj + teta_s*(Eigen::MatrixXd::Identity(b, b));
    Eigen::MatrixXd aux =  Rhj*Vj.inverse();

    Eigen::VectorXd aux_eij = aux*eij;
    u.row(j) = aux_eij;

    Eigen::MatrixXd sigma_aux_rhj = Sigmau - aux*t_Rhj;
    Eigen::VectorXd diag_var_u = sigma_aux_rhj.diagonal();
    var_u.row(j) = diag_var_u;

    Eigen::VectorXd zj_aux_eij = zj*aux_eij;

    Eigen::VectorXd yhat1 = xj_beta + zj_aux_eij;
    yhat.segment(a,b) = yhat1;

    Eigen::VectorXd vj = eij - zj_aux_eij;
    v.segment(a,b) = vj;

    //double auxil = teta_s*(1 - 1/b);
    //Eigen::VectorXd var2 = Eigen::VectorXd::Constant(b,auxil);
    //var_v.segment(a,b) = var2;
  }
  return Rcpp::List::create(Rcpp::Named("s_matc")=s_matc,
                            Rcpp::Named("s_matd")=s_matd,
                            Rcpp::Named("u")=u,
                            Rcpp::Named("var_u")=var_u,
                            Rcpp::Named("v")=v,
                            Rcpp::Named("yhat")=yhat);
}
