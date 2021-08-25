%NUMBER OF UNITS
N = 100;

%Define the connectivity matrix AKA ADJACENCY MATRIX
%Number of units

C = zeros(N);
for i=1:N-1
    C(i,i+1) = 1;
    C(i+1,i) = 1;
end
%% Define the distance matrix 
D = zeros(N);
x1 = 1;
for r=1:N
    x1 = x1 + 1;
    for c=x1:N
        D(r,c) = D(r,c-1) + 1;
        D(c,r) = D(r,c);
    end
end
%% Define parameters

% Node natural frequency in Hz 
          f=40;    % (i.e. f=40)):          

% Mean Delay in seconds (or speed)
          MD=0.012;  %(i.e. MD=0.01)             

% Global Coupling strength
          K=50; 
          
%% Call here the function of the Network Model

[Phases_Save, dt_save] = Kuramoto_Delays_Run(C,D,f,K,MD);

%%

lead_eigs = zeros(2,size(Phases_Save,2));
snd_eigs = zeros(2,size(Phases_Save,2));

for t=1:size(Phases_Save,2)
    
    % 2D data matrix 
    data_matrix = zeros(N,2);
    data_matrix(:,1) = cos(Phases_Save(:,t));
    data_matrix(:,2) = sin(Phases_Save(:,t));
    
    % return a matrix whose diagonals correspond to the covariance values
    % of each column 
    cov_mat = cov(data_matrix); 
    
    % EIGENDECOMPOSITION
    [vecs,vals] = eig(cov_mat);
    
    val1 = max(vals(:));
    val2 = second_max(vals); % -> Second highest eigenvalue
    
    % Determine the index of the highest eigenvector
    [~, col1] = find(vals==val1); 
    [~, col2] = find(vals==val2);
   
    %determine the eigenvectors
    vec1 = vecs(:,col1); 
    vec2 = vecs(:,col2); % -> Second highest eigenvalue
    
    %create matrices which store the leading and second eigenvectors of
    %each time instant
    lead_eigs(:,t) = lead_eigs(:,t) + vec1;
    snd_eigs(:,t) = snd_eigs(:,t) + vec2;
end

%% Define a new Order Parameter

OP = zeros(1,size(Phases_Save,2)); %row vector

for t=1:size(Phases_Save,2)
    
    % write the phasors coordinates
    x_coor = cos(Phases_Save(:,t).');
    y_coor = sin(Phases_Save(:,t).');
    phasors = [x_coor;y_coor];
    
    % take the dot product of all the phasors with the leading eig to get
    % the projections
    
    leading_eig = lead_eigs(:,t);
    products = zeros(1,size(Phases_Save,2));
    
    
    for col=1:N
        products(col) = dot(phasors(:,col),leading_eig);
    end
    
    OP(t) = mean(abs(products));
end
%% Video showing the first and second eigenvectors

  figure
  colormap jet
  
  B0=zeros(N,1);
  
  for t=1:size(Phases_Save,2)
      
      % All phasors
      subplot(2,2,1)
      cla
      hold on
      B1=Phases_Save(:,t);
      plot([B0 cos(B1)]',[B0 sin(B1)]')
      xlim([-1 1])
      ylim([-1 1])
      ylabel('Imag')
      xlabel('Real')
      axis square
      title('All phasors', {['t=' num2str(t*dt_save) 's']})
      
      % First and second leading eigenvectors
      subplot(2,2,2)
      cla
      hold on
      x1 = lead_eigs(1,t);
      y1 = lead_eigs(2,t);
      x2 = snd_eigs(1,t);
      y2 = snd_eigs(2,t);
      plot([0 x1]',[0 y1]',[0 x2]',[0 y2]')
      xlim([-1 1])
      ylim([-1 1])
      ylabel('Imag')
      xlabel('Real')
      axis square
      title('First and second leading eigenvectors')
      
      % Order parameter
      subplot(2,2,[3,4])
      plot(0:dt_save:(size(Phases_Save,2)-1)*dt_save,OP)
      % ylim([0 2])
      title('Order Parameter')
      pause(0.01)
      
  end
  
%% USEFUL FUNCTIONS

%Function to get the second highest eigenvalue
function[y] = second_max(x)
   y = max(x(x<max(x,[],'all')));
end