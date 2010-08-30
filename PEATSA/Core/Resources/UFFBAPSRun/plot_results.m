stability
#stability50

cstability = Astability\bstability;
eststability = Astability*cstability;

min_e = min(bstability);
max_e = max(bstability);

x = [min_e,0.1,max_e];

reg = [bstability ones(length(bstability),1)];
coeff = reg\eststability;

xlabel('Exp stability change');
ylabel('Calculated stability change');
plot(bstability, eststability, '*;data;',x,x,';optimal;',x,coeff(1)*x+coeff(2),';regression;');
print("stability.eps", "-depsc");




CORRCOEF_stability = corrcoef(bstability,eststability)
SPEARMAN_stability = spearman(bstability,eststability)
STANDARD_stability = std(eststability-bstability)
MEANUNSI_stability = sum(abs(eststability-bstability))/size(eststability-bstability)(1)
WANGSPMN_stability = 1-(6*sum(((ranks(bstability)-ranks(eststability)).^2)))/(rows(bstability)*(rows(bstability)^2-1))
