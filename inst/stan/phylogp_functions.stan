  // ---------------------------------------------------------------------------
  // phylogp: shared Stan functions.
  // Tree kernels of branch-wise form
  //   K = B diag(w) B' + rv 11' + diag(t2)
  // (w: edge variances, rv: root variance, t2: extra tip variance).
  // Edges are in ape "postorder" (children before parents), so a single pass
  // over the edges gives the O(N) pruning likelihood.
  // ---------------------------------------------------------------------------

  // Dense kernel matrix  K = B diag(w) B' + rv + diag(t2)
  matrix tree_kernel(matrix B, vector w, real rv, vector t2) {
    matrix[rows(B), rows(B)] K = diag_post_multiply(B, w) * B';
    return add_diag(K + rv, t2);
  }

  // O(N) Gaussian pruning log density of residuals r
  // (Felsenstein 1973; Ho & Ane 2014)
  real tree_pruning_lpdf(vector r, vector w, real rv, vector t2,
                         array[] int anc, array[] int des,
                         array[] int first_child, int n_nodes, int root) {
    int N = num_elements(r);
    int E = num_elements(w);
    vector[n_nodes] m;
    vector[n_nodes] v;
    real lp = 0;
    m[1:N] = r;
    v[1:N] = t2;
    for (e in 1:E) {
      int c = des[e];
      int p = anc[e];
      real vc = v[c] + w[e];
      if (first_child[e] == 1) {
        m[p] = m[c];
        v[p] = vc;
      } else {
        real s = v[p] + vc;
        lp += normal_lpdf(m[p] | m[c], sqrt(s));
        m[p] = (m[p] * vc + m[c] * v[p]) / s;
        v[p] = v[p] * vc / s;
      }
    }
    lp += normal_lpdf(0 | m[root], sqrt(v[root] + rv));
    return lp;
  }

  // Moments of the latent tip value x_i given y_{-i}, for every tip, by an
  // upward-downward pass of Gaussian messages (O(N)). Column 1: mean,
  // column 2: variance (including the tip edge, excluding t2[i]).
  // p(y_i | y_-i) = N(r_i | mean_i, var_i + t2_i).
  matrix tree_loo_moments(vector r, vector w, real rv, vector t2,
                          array[] int anc, array[] int des,
                          int n_nodes, int root) {
    int N = num_elements(r);
    int E = num_elements(w);
    vector[n_nodes] um = rep_vector(0, n_nodes);
    vector[n_nodes] uv = rep_vector(0, n_nodes);
    vector[n_nodes] P = rep_vector(0, n_nodes);
    vector[n_nodes] H = rep_vector(0, n_nodes);
    vector[n_nodes] om = rep_vector(0, n_nodes);
    vector[n_nodes] ov = rep_vector(0, n_nodes);
    matrix[N, 2] out;
    um[1:N] = r;
    uv[1:N] = t2;
    // upward: likelihood message from the subtree below each node
    for (e in 1:E) {
      int c = des[e];
      int p = anc[e];
      real vv;
      if (c > N) {
        uv[c] = 1 / P[c];
        um[c] = H[c] / P[c];
      }
      vv = uv[c] + w[e];
      P[p] += 1 / vv;
      H[p] += um[c] / vv;
    }
    // downward: message from everything outside the subtree (preorder)
    om[root] = 0;
    ov[root] = rv;
    for (k in 1:E) {
      int e = E - k + 1;
      int c = des[e];
      int p = anc[e];
      real vv = uv[c] + w[e];
      real m_ex;
      real v_ex;
      if (ov[p] <= 0) {          // parent value known exactly (fixed root)
        m_ex = om[p];
        v_ex = 0;
      } else {
        real Pp = P[p] - 1 / vv + 1 / ov[p];
        real Hp = H[p] - um[c] / vv + om[p] / ov[p];
        m_ex = Hp / Pp;
        v_ex = 1 / Pp;
      }
      om[c] = m_ex;
      ov[c] = v_ex + w[e];
    }
    out[, 1] = om[1:N];
    out[, 2] = ov[1:N];
    return out;
  }

  // Exact leave-one-out predictive log densities, dense Gaussian model
  // (Sundararajan & Keerthi 2001; Burkner, Gabry & Vehtari 2021)
  vector gaussian_loo_lpd(vector r, matrix K) {
    int N = num_elements(r);
    matrix[N, N] Kinv = inverse_spd(K);
    vector[N] g = Kinv * r;
    vector[N] out;
    for (i in 1:N) {
      real s2 = 1 / Kinv[i, i];
      out[i] = normal_lpdf(r[i] | r[i] - g[i] * s2, sqrt(s2));
    }
    return out;
  }

  // Mean and variance of r_i given r_{-i} for every i, dense Gaussian model
  // (the variance includes all diagonal terms of K)
  matrix gaussian_loo_moments(vector r, matrix K) {
    int N = num_elements(r);
    matrix[N, N] Kinv = inverse_spd(K);
    vector[N] g = Kinv * r;
    matrix[N, 2] out;
    for (i in 1:N) {
      real s2 = 1 / Kinv[i, i];
      out[i, 1] = r[i] - g[i] * s2;
      out[i, 2] = s2;
    }
    return out;
  }
