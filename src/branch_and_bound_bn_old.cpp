// // branch_and_bound_bn.cpp
// #include <Rcpp.h>
// #include <chrono>
// #include "branch_and_bound_bn.h"
// using namespace Rcpp;
//
// // Helper: convert integer to bitmask vector
// std::vector<int> intToBitsVec(int n, int p) {
//   std::vector<int> bits(p, 0);
//   for (int i = 0; i < p; i++) {
//     bits[i] = (n >> i) & 1;
//   }
//   return bits;
// }
//
// // [[Rcpp::export]]
// List branch_and_bound_bn_rcpp(int p,
//                              NumericMatrix X,
//                              List scores_list) {
//
//   int max_mask = 1 << p;
//   std::vector<int> global_best_parent_masks(p, 0);
//   double global_best_score = R_NegInf;
//   int nodes_explored = 0;
//   int pruned = 0;
//
//   NumericVector best_possible(p);
//   for (int j = 0; j < p; j++) {
//     NumericVector scores = scores_list[j];
//     double max_val = -R_PosInf;
//     for (int i = 0; i < scores.size(); i++) {
//       if (!NumericVector::is_na(scores[i]) && scores[i] > max_val)
//         max_val = scores[i];
//     }
//     best_possible[j] = max_val;
//   }
//
//   std::vector<int> prefix;
//   std::vector<int> parent_masks_so_far(p, 0);
//
//   std::function<void(std::vector<int>, std::vector<int>, double, int)> dfs;
//
//   dfs = [&](std::vector<int> prefix,
//             std::vector<int> parent_masks_so_far,
//             double current_score,
//             int remaining_nodes_mask) {
//
//     nodes_explored++;
//
//     if (remaining_nodes_mask == 0) {
//       if (current_score > global_best_score) {
//         global_best_score = current_score;
//         global_best_parent_masks = parent_masks_so_far;
//       }
//       return;
//     }
//
//     double optimistic = current_score;
//     for (int j = 0; j < p; j++) {
//       if ((remaining_nodes_mask >> j) & 1)
//         optimistic += best_possible[j];
//     }
//     if (optimistic <= global_best_score) {
//       pruned++;
//       return;
//     }
//
//     for (int j = 0; j < p; j++) {
//       if (!((remaining_nodes_mask >> j) & 1)) continue;
//
//       int preds_mask = 0;
//       for (int k : prefix) preds_mask |= (1 << k);
//
//       NumericVector scores_j = scores_list[j];
//       double best_local = -R_PosInf;
//       int chosen_mask = 0;
//
//       for (int sub = preds_mask;; sub = (sub - 1) & preds_mask) {
//         if (!NumericVector::is_na(scores_j[sub]) &&
//             scores_j[sub] > best_local) {
//           best_local = scores_j[sub];
//           chosen_mask = sub;
//         }
//         if (sub == 0) break;
//       }
//
//       if (best_local == -R_PosInf) continue;
//
//       std::vector<int> pm_new = parent_masks_so_far;
//       pm_new[j] = chosen_mask;
//       std::vector<int> prefix_new = prefix;
//       prefix_new.push_back(j);
//
//       int remaining_new = remaining_nodes_mask & (~(1 << j));
//
//       dfs(prefix_new, pm_new, current_score + best_local, remaining_new);
//     }
//   };
//
//   int all_nodes_mask = (1 << p) - 1;
//
//   // ============================
//   // ✔ FIXED TIMING USING std::chrono
//   // ============================
//   auto t0 = std::chrono::high_resolution_clock::now();
//
//   dfs(prefix, parent_masks_so_far, 0.0, all_nodes_mask);
//
//   auto t1 = std::chrono::high_resolution_clock::now();
//   double elapsed = std::chrono::duration<double>(t1 - t0).count();
//   // ============================
//
//   NumericMatrix best_adj(p, p);
//   for (int j = 0; j < p; j++) {
//     int mask = global_best_parent_masks[j];
//     for (int i = 0; i < p; i++) {
//       if ((mask >> i) & 1)
//         best_adj(i, j) = 1;
//     }
//   }
//
//   return List::create(
//     Named("best_adj") = best_adj,
//     Named("best_score") = global_best_score,
//     Named("time") = elapsed,
//     Named("nodes_explored") = nodes_explored,
//     Named("pruned") = pruned
//   );
// }
//
//
