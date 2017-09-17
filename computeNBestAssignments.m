%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computes the N most likely assignments given an assignment costMatrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function assignList = computeNBestAssignments(n, costMatrix, costUnassignedTracks, costUnassignedDetections)
% Compute the N best assignments, by modifying the cost matrix after
% finding the N-1 best assignments. 
% This would be so much easier/faster if non-assignments were not possible :(

% assignment is a struct: matches, unassignedUAVs, unassignedDetections

assignList = struct('matches',repmat({nan}, 1,n), 'unassignedUAVs',cell(1,n), 'unassignedDetections',cell(1,n), 'score',cell(1,n));
nAgents = size(costMatrix, 1);
assignMatrix = nan(nAgents, n);

for i = 1:n
    if all(isnan(assignMatrix))
        newAssignment = struct;
        [m, a, b] = assignDetectionsToTracks(-log(costMatrix), costUnassignedTracks, costUnassignedDetections);
        newAssignment.matches = m;
        newAssignment.unassignedUAVs = a;
        newAssignment.unassignedDetections = b;
        newAssignment.score = computeAssignmentScore(newAssignment, costMatrix, costUnassignedTracks, costUnassignedDetections);
    else
        newAssignment = computeBestAssignmentExcluding(nAgents, costMatrix, costUnassignedTracks, costUnassignedDetections, assignMatrix);
    end
    % update the prevAssignments
    assignList(i) = newAssignment;
 
    sortedAssignment = sortrows([newAssignment.matches; [newAssignment.unassignedUAVs(:) NaN(length(newAssignment.unassignedUAVs),1)]]);
    assignMatrix(:,i) = sortedAssignment(:, 2);
    %for assignRow = newAssignment.matches',
    %    assignMatrix(assignRow(2), i) = assignRow(1); 
    %end        
end

end


function [score] = computeAssignmentScore(assignInfo, costMatrix, costUnassignedTracks, costUnassignedDetections)
    nUnassigned = length(assignInfo.unassignedUAVs);
    nUndetected = length(assignInfo.unassignedDetections);
    score = prod(costMatrix(sub2ind(size(costMatrix), assignInfo.matches(:,1), assignInfo.matches(:,2))))...
        *(costUnassignedTracks^-nUnassigned)*(costUnassignedDetections^-nUndetected);

end


function [bestAssignment] = computeBestAssignmentExcluding(P, prevCostMatrix, costUnassignedTracks, costUnassignedDetections, prevAssignments)
% this requires some dynamic programming to make it really efficient: we
% need to store previously considered subspace scores, as they may be
% better than the new scores we generate

% prevAssignments is an PxQ matrix, where P is the number of agents to
% assign, and Q is the number of previous assignments made.
% prevAssignments(P,A) = the assignment given to agent P in round A
% If an agent was unassigned, the value will be nan - should still work

% savedAssignments is a struct array of the N-1(?) best scoring assignments 
% from the previous iteration, including scores

 bestScore = -inf; % best assignment so far - look in savedAssignments
 bestAssignment = nan;
 
 for i = 1:P
     % for each agent, we exclude its previous assignments by making those
     % entries in the costMatrix humunguous
     
     agentAssignments = prevAssignments(i,:);
     currentCostMatrix = prevCostMatrix;
     goodIdxs = agentAssignments(~isnan(agentAssignments));
     currentCostMatrix(i, goodIdxs) = 1e-28;
     %for j = goodIdxs,
     %    currentCostMatrix(i,j) = currentCostMatrix(i,j)*1e-5;
     %end


     %[matches, unUAV, unDet] = assignDetectionsToTracks(-log(currentCostMatrix), costUnassignedTracks, costUnassignedDetections);
     [matches, unUAV, unDet] = assignDetectionsToTracks(-log(currentCostMatrix), 50000, 50000);
     assignment.matches = matches; assignment.unassignedUAVs = unUAV; assignment.unassignedDetections = unDet; 
     assignment.score = computeAssignmentScore(assignment, prevCostMatrix, costUnassignedTracks, costUnassignedDetections); % see how good this assignment is
     if assignment.score > bestScore
         % so far, this is the best we've got
         bestAssignment = assignment;
         bestScore = assignment.score;
     end
     
     
 end
end
             