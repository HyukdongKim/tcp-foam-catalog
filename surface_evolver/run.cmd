// Surface Evolver protocol: relax, refine (r), relax ... ; report cost = 2*area/N at each level
procedure report(integer lev) {
  printf "LEVEL %d facets %d cost %.10f maxVolErr %.1e\n", lev, facet_count, 2*total_area/body_count, max(body, abs(volume - target));
};
procedure settle() {
  local prev, it;
  g 100; u; V; g 100;
  prev := 1e30; it := 0;
  while (abs(prev - total_area) > 1e-10*total_area) and (it < 80) do {
    prev := total_area; g 100; u; hessian_seek; it := it + 1;
  };
  printf "  settle iterations %d\n", it;
};
settle(); report(0);
r; settle(); report(1);
r; settle(); report(2);
r; settle(); report(3);
q
