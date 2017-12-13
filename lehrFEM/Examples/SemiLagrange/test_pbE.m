function test_pbE(NewMesh,pMesh,A);

nEdges=size(NewMesh.Edges,1);

for i=1:nEdges
    ip=[sum(A(i,:)'.*(NewMesh.Coordinates(NewMesh.Edges(:,2),1)-NewMesh.Coordinates(NewMesh.Edges(:,1),1))) ...
        sum(A(i,:)'.*(NewMesh.Coordinates(NewMesh.Edges(:,2),2)-NewMesh.Coordinates(NewMesh.Edges(:,1),2)))];
    e=pMesh.Coordinates(pMesh.Edges(i,2),:)-pMesh.Coordinates(pMesh.Edges(i,1),:);
    %norm(ip-e)
    if norm(ip-e)>=eps
        [i ip e]
    end
    
end