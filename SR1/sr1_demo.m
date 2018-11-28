disp('Loading example data');
load ex.mat;

Res = DATA;

disp('Calculating shifted rank-1 approximation using 3 matrices');

for k=1:3
    disp(['... (' num2str(k) '/3)']);
    
    % calculate approximation of resudal
    [B{k},u{k},v{k},lambda{k}] = sr1(Res);
    
    % update residual
    Res = Res-B{k};
    
    % for better visuability shift each event such that its peak is at k/4
    % of its length
    [~,pos] = max(abs(u{k}));
    u{k} = circshift(u{k},round(pos+size(DATA,2)/4*k));
    
    % also shift the vector lambda to show the actual position of the peak.
    % We have do calculate position modulo data length since the method
    % uses circulant shifts
    lambda{k} =  mod(lambda{k} + pos,size(DATA,2));
end

disp('Plotting results');
figure;

% FIRST SUBPLOT: original data with shift vectors indicating the peak of
% all events
subplot(231); imagesc(DATA); colormap pink; hold on;
plot(lambda{1},'x'); plot(lambda{2},'x'); plot(lambda{3},'x'); title('Data with tracked events');

% SECOND SUBPLOT: extracted waveforms
subplot(232); hold on; plot(u{1}); plot(u{2}); plot(u{3}); axis tight; legend({'Event 1','Event 2','Event 3'}); title('Extracted waveforms');

% THIRD SUBPLOT: amplitudes of all events
subplot(233); hold on; plot(abs(v{1})); plot(abs(v{2})); plot(abs(v{3})); axis tight; legend({'Event 1','Event 2','Event 3'}); title('Event amplitudes');

% SUBPLOT 4-6: all events
subplot(234); imagesc(B{1}); colormap pink; title('Event 1');
subplot(235); imagesc(B{2}); colormap pink; title('Event 2');
subplot(236); imagesc(B{3}); colormap pink; title('Event 3');

disp('This example shows two drawbacks of the method that will be considered in future work:');
disp('- Overlapping events are hard to seperate and cause disturbance in the data');
disp('- Events that locally fanish (e.g., Event 3) are not yet implemented. Use e.g., a threshold on the amplitude to account these phenomenom');