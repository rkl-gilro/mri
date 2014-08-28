function plot_triangulations(tr1, tr2)
%% Plot triangulations %%

%% Plot the original triangulation and reference points.
figure;
subplot(1,2,1);
tetramesh(tr1); title('First mesh'); hold on;
alpha(.1)

subplot(1,2,2);
tetramesh(tr2);title('Second mesh'); 
alpha(.1)
hold on;

axis equal;