function threeDPrinter(action, xPos, zPos)
% THREEDPRINTER  A simple controller for a modified 3D printer (X and Z only).
%                Uses a persistent serial port connection.
%
% Usage:
%   threeDPrinter("open");
%       - Open the serial port; optionally set up G-code mode (e.g., G90).
%
%   threeDPrinter("home");
%       - Home the X and Z axes (G28 X Z).
%
%   threeDPrinter("move", x, z);
%       - Move the printer to (x, z) in absolute coordinates.
%
%   threeDPrinter("close");
%       - Closes the serial port.
%
% Example in a script:
%   threeDPrinter("open");
%   threeDPrinter("home");
%
%   for i = 1:10
%       threeDPrinter("move", 10*i, 5);   % move in X, fix Z=5
%       % capture radar data ...
%   end
%
%   threeDPrinter("close");

    persistent s  % The serial port handle persists across calls

    % ------------------- ADJUST THESE TO YOUR SETUP -------------------
    comPort  = "COM5";    % e.g. "COM3" on Windows, or "/dev/ttyUSB0" on Linux
    baudRate = 115200;    % or 250000, etc.
    feedRate = 3000;      % mm/min feed rate for the G1 moves
    % -----------------------------------------------------------------

    % Decide which action to take
    switch lower(action)
        %===============================================================
        case "open"
            % Only open if not already open
            if ~isempty(s) && isvalid(s)
                disp("[threeDPrinter] Serial port already open.");
                return;
            end

            disp("[threeDPrinter] Opening serial port...");
            s = serialport(comPort, baudRate);

            % Small pause to let the printer/firmware settle
            pause(2);

            % Flush any initial data in the buffer
            flush(s);

            % Switch to absolute positioning by default
            writeline(s, "G90");  
            pause(0.5);

            disp("[threeDPrinter] Port open and ready.");

        %===============================================================
        case "home"
            % Ensure the port is open
            if isempty(s) || ~isvalid(s)
                error("[threeDPrinter] Port not open. Call threeDPrinter('open') first.");
            end

            disp("[threeDPrinter] Homing X and Z axes...");
            writeline(s, "G28 X Z");  % Home X and Z
            % You may need a longer pause depending on your machine
            pause(5);  
            disp("[threeDPrinter] Homing complete.");

        %===============================================================
        case "move"
            % Ensure the port is open
            if isempty(s) || ~isvalid(s)
                error("[threeDPrinter] Port not open. Call threeDPrinter('open') first.");
            end

            % Validate inputs
            if nargin < 3
                error('[threeDPrinter] "move" requires xPos, zPos inputs.');
            end

            % Build G-code command
            gcodeMove = sprintf("G1 X%.2f Z%.2f F%d", xPos, zPos, feedRate);
            writeline(s, gcodeMove);

            % Short pause to ensure the printer starts moving;
            % Adjust if you need guaranteed completion before next command.
            pause(0.5);

        %===============================================================
        case "close"
            % Close (clear) the port if open
            if ~isempty(s)
                disp("[threeDPrinter] Closing serial port...");
                clear s;  % or delete(s);
                s = [];
            else
                disp("[threeDPrinter] Port is not open.");
            end

        %===============================================================
        otherwise
            error('[threeDPrinter] Invalid action "%s". Use "open", "home", "move", or "close".', action);
    end
end
