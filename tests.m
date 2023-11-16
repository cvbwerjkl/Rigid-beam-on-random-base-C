%based on A.J.M. Ferreira book MATLAB Codes for FEA

printf("Octave start\n");
deviation = 0.15;
lBeam = 25;
corRadius = 1;
EI = 36 * 100 * 40 * 40 * 40 / 12;
P = -5;
numberSprings = 2; 
springsVector = [1000 2000];

numberElements = numberSprings * 2;
nodeCoordinates = linspace(0,lBeam,numberElements+1)';
numberNodes = size(nodeCoordinates,1);
xx = nodeCoordinates(:,1);
elementNodes = zeros(numberElements,2);
for i = 1:numberElements
    elementNodes(i,1) = i; 
    elementNodes(i,2) = i+1;
end

GDof = 2*numberNodes; 

force = zeros(GDof,1);
stiffness = zeros(GDof);
% calculation of the system stiffness matrix
% and force vector
for e = 1:numberElements
    % elementDof: element degrees of freedom (Dof)
    indice = elementNodes(e,:);
    elementDof = [2*(indice(1)-1)+1 2*(indice(2)-1) ...
        2*(indice(2)-1)+1 2*(indice(2)-1)+2];
    % length of element
    LElem = xx(indice(2))-xx(indice(1));
    k1 = EI/(LElem)^3*[12   6*LElem -12 6*LElem;
        6*LElem 4*LElem^2 -6*LElem 2*LElem^2;
        -12 -6*LElem 12 -6*LElem ;
        6*LElem 2*LElem^2 -6*LElem 4*LElem^2];
    
    f1 = [P*LElem/2 P*LElem*LElem/12 P*LElem/2 ...
        -P*LElem*LElem/12]';
    
    % equivalent force vector
    force(elementDof) = force(elementDof) + f1;
    
    % stiffness matrix
    stiffness(elementDof,elementDof) = ...
        stiffness(elementDof,elementDof)+k1;
end

% adding springs to stiffness matrix
springsMatrix = stiffness;
springsMatrix(3,3) += springsVector(1);
springsMatrix(7,7) += springsVector(2);

solution = springsMatrix\force;

% save forceM_Octave.txt force;
% save stiffnessM_Octave.txt stiffness;
% save stiffnessMS_Octave.txt springsMatrix;
save solution_Octave.txt solution;

printf("Octave done\n");