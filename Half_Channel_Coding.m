clear all
close all
clc

% to read the video
obj = VideoReader('highway.avi');
a = read(obj);

% to get the number of frames
frames = get(obj,'NumFrames');
trellis_half = poly2trellis(7,[171 133]);

% to extract the frames of the video to work on them
for i = 1 : frames
    I(i).cdata = a(:,:,:,i);
end


% to generate a new video with the same size as the original video
s = size(I(1).cdata);
mov(1:frames) = struct('cdata',zeros(s(1),s(2),3,'uint8'),'colormap',[]);

%probabilities to work with
p001 = 0.001;
p1 = 0.1;

% to extract the data of each colour in each frame
for i = 1 : frames
    % Red Components of the Frame
    R(:,:,i) = I(i).cdata(:,:,1);
    
    % Green Components of the Frame
    G(:,:,i) = I(i).cdata(:,:,2);
    
    % Blue Components of the Frame
    B(:,:,i) = I(i).cdata(:,:,3);
end

%     to convert the data from unsigned integers (the original format) to binary
%     takes two conversions to convert
    
    % from unsigned integers to double 
    Rdouble = double(R);
    Gdouble = double(G);
    Bdouble = double(B);
    
    % from double to binary
    Rbin = de2bi(Rdouble);
    Gbin = de2bi(Gdouble);
    Bbin = de2bi(Bdouble);
    
    %reshape frame into 1D array
    RbinR = reshape(Rbin,1,202752*30); % 144*176*8 = 202,752
    GbinR = reshape(Gbin,1,144*176*8*30);
    BbinR = reshape(Bbin,1,144*176*8*30);

    % concatenating all the colours together to make the full data
    fullFrame = [RbinR GbinR BbinR];
    
    % takes one dimentional array so feed fullFrame then split
    encoded = convenc(fullFrame, trellis_half);

    % passing the encoded through the channel
    error_half = bsc(encoded,p1);

    % decoding the data
    decoded = vitdec(error_half,trellis_half,35,'trunc','hard');
    
    % splitting the decoded data into packets
    decoded_packets = reshape(decoded,17820,1024); %rows = 17820, packet size = 1024
    
    % putting the packets together
    rec_frames = reshape(decoded_packets,1,144*176*8*3*30);

    % splitting the frames by colour
    recRbin = decoded_packets(1:6082560);
    recGbin = decoded_packets(6082560+1:6082560*2);
    recBbin = decoded_packets(6082560*2+1:6082560*3);

    % reshaping from a 1-D array to a 2-D as 144x176 by 8
    recRbinR = reshape(recRbin,144*176*30,8);
    recGbinR = reshape(recGbin,144*176*30,8);
    recBbinR = reshape(recBbin,144*176*30,8);
    
    % binary to double
    recRdbl = bi2de(recRbinR);
    recGdbl = bi2de(recGbinR);
    recBdbl = bi2de(recBbinR);
    
    % double to unsigned integer 8
    recRint = uint8(recRdbl);
    recGint = uint8(recGdbl);
    recBint = uint8(recBdbl);
    
    % reshaping 144x176 flat to 176 by 144 times 30 for the 30 frames
    recRre = reshape(recRint,144,176,30);
    recGre = reshape(recGint,144,176,30);
    recBre = reshape(recBint,144,176,30);
    
    % unloading the frames into the empty video created before
    for i = 1 : frames
        mov(1,i).cdata(:,:,1) = recRre(:,:,i);
        mov(1,i).cdata(:,:,2) = recGre(:,:,i);
        mov(1,i).cdata(:,:,3) = recBre(:,:,i);
    end

    % creating the video file
    obj2 = VideoWriter('after_half_p1.avi')
    open(obj2)
    writeVideo(obj2,mov)
    close(obj2)

    implay('after_half_p1.avi')
    
    % probabilities for testing
    probz = [0.0001:0.02:0.2];
    
    for i = 1 : length(probz)
        %takes one dimentional array so feed fullFrame then split
        encoded = convenc(fullFrame, trellis_half);
        error_half = bsc(encoded,probz(i));
        decoded = vitdec(error_half,trellis_half,35,'trunc','hard');
        decoded_packets = reshape(decoded,17820,1024);
        
        % gathering the BER
        [number(i),BER(i)] = biterr(fullFrame,decoded)
    end
    
    %plotting the values
    figure()
    plot(probz,BER)
    xlabel("Probability of Error")
    ylabel("Bit Error Rate")
    title("Half Channel Coding BER")

