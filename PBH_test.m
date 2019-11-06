close all;
%%
lmda = eig(Asys);
rnk_pbh = zeros(size(lmda));
for i = 1:length(lmda)
    rnk_pbh(i) = rank([Asys - lmda(i)*eye(2*n), Bcont],1e-5);
end

[Va, Da] = eig(Asys);

%% Gramian
sys = ss(Asys,Bcont,Csys,0);
WC = gram(sys,'c');
WO = gram(sys,'o');
[Ua,Sa,~] = svd(WC);
Sa_diag = diag(Sa);
%% Plots
for i=1:5
    figure;
    plot(real(Va(:,i)/max(Va(:,i))))
    hold on;
    plot(Ua(:,i)/max(Ua(:,i)));hold off;
end