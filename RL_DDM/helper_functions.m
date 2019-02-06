% Put the results of the model (learning rate) into the data struct

for i=1:length(A)
    A(i).lr = A_model.x(i,4)
end 

% Making a large list of reaction times paired with kmeans learning rate
% assignments.

for i=1:length(A)
    for j=1:length(A(i).rt)
        rt_list(i,j,1) = A(i).rt(j);
        rt_list(i,j,2) = A(i).lr;
        rt_list(i,j,3) = A(i).kk;
    end
end

    
