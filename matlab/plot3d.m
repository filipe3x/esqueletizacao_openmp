[W,P] = meshgrid(0:256:8192,1:1:32);
Z = (W./sqrt (2)) .* ( 2.*(3 + W./500) + 2.*3.*(log (P)) ) + ((W.*W) ./ P)./500 + 10 + (((W.*W) ./ P)./500  + 10) .* (P - 1);

%# plot 3D bars
A=axes;
b = surf(W,P,Z);
colormap(jet);
colorbar
xlabel('Width'), ylabel('Procs'), zlabel('microseconds')
%set(gca, 'XTickLabel',W, 'YTickLabel',P)
%set(A,'XScale','log');
%set(A,'YScale','log');

