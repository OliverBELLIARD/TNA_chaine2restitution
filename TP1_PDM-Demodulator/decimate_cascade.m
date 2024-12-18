function Xout = decimate_cascade(Xin, Fs, M1, M2, M3, M4, M5, M6, M7)
%% Decimates in a cascade of up to 7 stages

%% Cascade Decimation
% Stage 1 - Decimate by M1
Fs1 = Fs / M1;
pdm_filtered1 = decimate_audio(Xin, Fs, M1); % Pass original Fs to Stage 1

if nargin > 3
    % Stage 2 - Decimate by M2
    Fs2 = Fs1 / M2;
    pdm_filtered2 = decimate_audio(pdm_filtered1, Fs1, M2); % Pass Fs1 to Stage 2

    if nargin > 4
        % Stage 3 - Decimate by M3
        Fs3 = Fs2 / M3;
        pdm_filtered3 = decimate_audio(pdm_filtered2, Fs2, M3); % Pass Fs2 to Stage 3

        if nargin > 5
            % Stage 4 - Decimate by M4
            Fs4 = Fs3 / M4;
            pdm_filtered4 = decimate_audio(pdm_filtered3, Fs3, M4); % Pass Fs3 to Stage 4

            if nargin > 6
                % Stage 4 - Decimate by M4
                Fs5 = Fs4 / M5;
                pdm_filtered5 = decimate_audio(pdm_filtered4, Fs4, M5); % Pass Fs3 to Stage 4

                if nargin > 7
                    % Stage 6 - Decimate by M6
                    Fs6 = Fs5 / M6;
                    pdm_filtered6 = decimate_audio(pdm_filtered5, Fs5, M6); % Pass Fs3 to Stage 4

                    if nargin > 8
                        % Stage 7 - Decimate by M7
                        %Fs7 = Fs6 / M7;
                        Xout = decimate_audio(pdm_filtered6, Fs6, M7); % Pass Fs3 to Stage 4
                    else
                        Xout = pdm_filtered6;
                    end
                else
                    Xout = pdm_filtered5;
                end
            else
                Xout = pdm_filtered4;
            end
        else
            Xout = pdm_filtered3;
        end
    else
        Xout = pdm_filtered2;
    end
else
    Xout = pdm_filtered1;
end

end