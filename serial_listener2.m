% serial_listener.m
port = "COM3";
baudRate = 115200;
s = serialport(port, baudRate);
configureTerminator(s, "CR/LF");   % Arduino println uses CR/LF

disp("Listening...");
while true
    if s.NumBytesAvailable > 0
        line = strtrim(readline(s));

        % Debug: show what MATLAB received
        fprintf("Received from Arduino: %s\n", line);

        if startsWith(line, "RUN:")
            parts = split(extractAfter(line,"RUN:"), ",");
            D  = str2double(parts(1));
            d  = str2double(parts(2));
            U0 = str2double(parts(3));

            % Debug: show parsed values
            fprintf("Parsed D = %.0f mm, d = %.0f mm, U0 = %.2f m/s\n", D, d, U0);

            % Run your simulation
            [P1, P2] = simulation(D, d, U0);

            % Debug: show what MATLAB is sending back
            fprintf("Sending to Arduino: %.2f, %.2f\n", P1, P2);

            writeline(s, sprintf('%.2f,%.2f', P1, P2));
        end
    end
    pause(0.01);
end
