%Nested for loop to create grid for 3D Printer
%Nested for loop is better because vertical speed is slower than horizontal
%move in zigzag pattern to eliminate unnecessary movement of row reset
%pause line included to account for varying measurement ties
%9 second pause for vertical movement
%24 second pause from home to starting position

Length = 200; % mm
Height = 200; % mm
step_size = 40; % mm

% Initialize
threeDPrinter("open");
threeDPrinter("home");
threeDPrinter("move", 50, 100);
pause(24);

% Manually move through the 5x5 grid
threeDPrinter("move", 90, 100);
pause(1);
threeDPrinter("move", 130, 100);
pause(1);
threeDPrinter("move", 170, 100);
pause(1);
threeDPrinter("move", 210, 100);
pause(1);
threeDPrinter("move", 250, 100);
pause(1);
threeDPrinter("move", 250, 140);
pause(9);
threeDPrinter("move", 210, 140);
pause(1);
threeDPrinter("move", 170, 140);
pause(1);
threeDPrinter("move", 130, 140);
pause(1);
threeDPrinter("move", 90, 140);
pause(1);
threeDPrinter("move", 50, 140);
pause(1);
threeDPrinter("move", 50, 180);
pause(9);
threeDPrinter("move", 90, 180);
pause(1);
threeDPrinter("move", 130, 180);
pause(1);
threeDPrinter("move", 170, 180);
pause(1);
threeDPrinter("move", 210, 180);
pause(1);
threeDPrinter("move", 250, 180);
pause(1);
threeDPrinter("move", 250, 220);
pause(9);
threeDPrinter("move", 210, 220);
pause(1);
threeDPrinter("move", 170, 220);
pause(1);
threeDPrinter("move", 130, 220);
pause(1);
threeDPrinter("move", 90, 220);
pause(1);
threeDPrinter("move", 50, 220);
pause(1);
threeDPrinter("move", 50, 260);
pause(9);
threeDPrinter("move", 90, 260);
pause(1);
threeDPrinter("move", 130, 260);
pause(1);
threeDPrinter("move", 170, 260);
pause(1);
threeDPrinter("move", 210, 260);
pause(1);
threeDPrinter("move", 250, 260);
pause(1);

threeDPrinter("home");

threeDPrinter("close");
