%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------
% Created by isardSAT S.L.
% ------------------------------------------------------
% Pablo N García
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Tibet_Lakes=kml2lla(kmlfile)

%kml_in=importdata('Qinghai Tibet Lakes.kml', ',');

kml_in=importdata(kmlfile);


%%%%%%%%   names   %%%%%%%%

% finding names
recsfound=0; 
for i=1:size(kml_in,1)
    if strfind(char(kml_in(i,1)),'name')
        recsfound=recsfound+1; Tibet_Lakes(recsfound,1).name=char(kml_in(i,1)); 
    end
end

% cleaning names
for i=1:recsfound
    string_1(i,:)=strfind(Tibet_Lakes(i,1).name,'<');
    string_2(i,:)=strfind(Tibet_Lakes(i,1).name,'>');
    Tibet_Lakes(i,1).name(string_1(i,2):string_2(i,2))=[];
    Tibet_Lakes(i,1).name(1:string_2(i,1))=[];
end

% deleting non useful names
Tibet_Lakes(1:2,:)=[];



%%%%%%%%   coordinates data   %%%%%%%%

%% finding coordinates
recsfound=0;
for i=1:size(kml_in,1)
    if strfind(char(kml_in(i,1)),'<coordinates>')
        recsfound=recsfound+1;
        temp_data=char(kml_in(i+1,1));
        initial_spaces = strfind(temp_data,'	');
        getcoord(recsfound,1).data=temp_data(max(initial_spaces)+1:size(temp_data,2)); %delete first 7 tabs
        getcoord(recsfound,1).data=changem(getcoord(recsfound,1).data,',0,',',0 '); % substitute ',0 ' by ',0,'
    end
end

%% re-ordering by lon / lat / alt

% find commas positions
for i=1:recsfound
    aa=char(getcoord(i,1).data(1,:));
    commas(1,i).pos=strfind(char(getcoord(i,1).data(1,:)),','); 
end

% separating lat, lon & alt data
for i=1:recsfound
    dataread=char(getcoord(i,1).data);
    data(i).coord(1,1)=str2num(dataread(1,1:(commas(1,i).pos(1,1)-1))); % first coordinate is before a comma
    for i_comma=1:(size(commas(1,i).pos,2)-1)
        data(i).coord(floor((i_comma)/3)+1,mod(i_comma,3)+1)=str2num(dataread(1,(commas(1,i).pos(1,i_comma)+1):(commas(1,i).pos(1,i_comma+1)-1))); % get every record between commas
    end
end

for i=1:recsfound
    Tibet_Lakes(i,1).coord=data(i).coord;
end





