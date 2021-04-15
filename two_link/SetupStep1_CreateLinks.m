RelativeBase = [0; 0; 0];        
RelativeFollower = [0; 0; 0.6];   
RelativeCoM = [0; 0; 0.3]; 
Mass = 20;
Inertia = eye(3);
Name = 'Shin';
save('datafile_Shin', 'RelativeBase', 'RelativeFollower', 'RelativeCoM', 'Mass', 'Inertia', 'Name');

RelativeBase = [0; 0; 0];        
RelativeFollower = [0; 0; 0.6];   
RelativeCoM = [0; 0; 0.3];  
Mass = 20;
Inertia = eye(3);
Name = 'Hip';
save('datafile_Hip', 'RelativeBase', 'RelativeFollower', 'RelativeCoM', 'Mass', 'Inertia', 'Name');

RelativeBase = [0; 0; 0];        
RelativeFollower = [0; 0; 0.6];   
RelativeCoM = [0; 0; 0.3]; 
Mass = 20;
Inertia = eye(3);
Name = 'Torso';
save('datafile_Torso', 'RelativeBase', 'RelativeFollower', 'RelativeCoM', 'Mass', 'Inertia', 'Name');


