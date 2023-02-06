clear all
close all
clc

% to read the video
obj = VideoReader('highway.avi');
a = read(obj);

% to get the number of frames
frames = get(obj,'NumFrames');

trellis = poly2trellis(7,[171 133]);

% code rate labels for tracing
code_rate_label=["8/9"	"4/5" "2/3" "4/7" "1/2"];

% code rates
code_rate = [1 1 1 0 1 0 1 0 0 1 1 0 1 0 1 0; 
			1 1 1 0 1 0 1 0 1 1 1 0 1 0 1 0;
			1 1 1 0 1 1 1 0 1 1 1 0 1 1 1 0;
			1 1 1 1 1 1 1 0 1 1 1 1 1 1 1 0;
			1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];

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
disp("after p")

% to extract the data of each colour in each frame
for i = 1 : frames
    % Red Components of the Frame
    R(:,:,i) = I(i).cdata(:,:,1);
    
    % Green Components of the Frame
    G(:,:,i) = I(i).cdata(:,:,2);
    
    % Blue Components of the Frame
    B(:,:,i) = I(i).cdata(:,:,3);
end

%   to convert the data from unsigned integers (the original format) to binary
%   takes two conversions to convert
    
    % from unsigned integers to double    
    Rdouble = double(R);
    Gdouble = double(G);
    Bdouble = double(B);
    
    % from double to binary
    Rbin = de2bi(Rdouble);
    Gbin = de2bi(Gdouble);
    Bbin = de2bi(Bdouble);
    
    % reshape frame into 1D array
    RbinR = reshape(Rbin,1,202752*30); % 144*176*8 = 202,752 * 30 frames
    GbinR = reshape(Gbin,1,144*176*8*30);
    BbinR = reshape(Bbin,1,144*176*8*30);

    % concatenating all the colours together to make the full data
    fullFrame = [RbinR GbinR BbinR];

    % code rate counter
    j = 1;

    % Encoding & Decoding the first time for the checking of the while loop %

    % encoding at a rate of 1/2
    encoded_half = convenc(fullFrame, trellis,code_rate(5,:));

    % encoding at a rate of 8/9
    % takes one dimentional array so feed fullFrame then split
    encoded = convenc(encoded_half, trellis,code_rate(j,:)); 

    % passing the encoded through the channel
  	error = bsc(encoded,p001);

    % decoding the data at rate 8/9
    decoded_half = vitdec(error,trellis,35,'trunc','hard',code_rate(j,:));
    
    % decoding the data at rate 1/2
    decoded = vitdec(decoded_half,trellis,35,'trunc','hard',code_rate(5,:));
    
    % splitting the decoded data into packets
    % rows = 17820, packet = 1024
 	decoded_packets = reshape(decoded,17820,1024); 

    % xor fullFrame and decoded to get the bits in error in an matrix
    % reshaping to flatten the array into 1-D
    % nnz checks the perimeter of a matrix so we give it a 1-D array
    % to check on the number of 1's
    % if the number of ones is greater than 0 then bits are in error
    % so we must change the code rate
    while nnz(reshape(xor(fullFrame,decoded).',1,[])) > 0
        encoded_half = convenc(fullFrame, trellis,code_rate(5,:));
        encoded = convenc(encoded_half, trellis,code_rate(j,:)); %takes one dimentional
        %array so feed fullFrame then split
        error = bsc(encoded,p001);
        decoded_half = vitdec(error,trellis,35,'trunc','hard',code_rate(j,:));
        decoded = vitdec(decoded_half,trellis,35,'trunc','hard',code_rate(5,:));
        
        % breaks out of loop if code rate reaches 1/2
      	if j == length(code_rate(:,1))
            break;
        else
          	disp("The following code rate has failed : " + code_rate_label(j))
      		j = j + 1;
        end
    end
     
    disp("The final code rate is : " + code_rate_label(j))
    
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
    obj2 = VideoWriter('after_incr_p001.avi')
    open(obj2)
    writeVideo(obj2,mov)
    close(obj2)

    implay('after_incr_p001.avi')

    %________________________________________________________

    % probabilities for obtaing the Bit Error Rate

    probz = [0.0001:0.04:0.2];
    
    % resetting the j value to 1
    j = 1;

    % Encoding & Decoding the first time for the checking of the while loop %
    
    % encoding at a rate of 1/2
    encoded_half = convenc(fullFrame, trellis,code_rate(5,:));

    % encoding at a rate of 8/9
    % takes one dimentional array so feed fullFrame then split
    encoded = convenc(encoded_half, trellis,code_rate(j,:));

    % passing the encoded through the channel
  	error = bsc(encoded,p001);

    % decoding the data at rate 8/9
    decoded_half = vitdec(error,trellis,35,'trunc','hard',code_rate(j,:));

    % naming it decoded1 to reset the check when continuing in the loop
    decoded1 = vitdec(decoded_half,trellis,35,'trunc','hard',code_rate(5,:));
    
    %rows = 17820, packet  = 1024 bits
 	decoded_packets = reshape(decoded1,17820,1024);
    disp("after decoding the first time")

    i = 1
    code_rate_size = length(code_rate(:,1))
    check = nnz(reshape(xor(fullFrame,decoded1).',1,[]))

    for i = 1 : length(probz)
        disp("iteration = " + i)
        while check > 0
            encoded_half = convenc(fullFrame, trellis,code_rate(5,:));
            encoded = convenc(encoded_half, trellis,code_rate(j,:)); %takes one dimentional array so feed fullFrame then split
            error = bsc(encoded,probz(i)); 

            decoded_half = vitdec(error,trellis,35,'trunc','hard',code_rate(j,:));
            decoded = vitdec(decoded_half,trellis,35,'trunc','hard',code_rate(5,:));
            disp("loop" + j)
            
            % updating the check value
            check = nnz(reshape(xor(fullFrame,decoded).',1,[]))
            
            if j == length(code_rate(:,1))
                j = 1
                % saving the BER values
                [ER(i),BER(i)] = biterr(fullFrame,decoded)
%                 i = i + 1
                continue
                
            else
                disp("The following code rate has failed : " + code_rate_label(j))
                disp("old j = " + j)
                j = j + 1
            end
        end
        
        % updating the check value
        check = nnz(reshape(xor(fullFrame,decoded1).',1,[]))
        disp(i)
    end
     
    disp("BER Completed")
    %-----------------------

    % generating the Throughput values
    for i = 1 : length(BER)
        Throughput(i) = 1 - BER(i)
    end
    
    %plotting the BER values
    figure()
    plot(probz,BER)
    xlabel("Probability of Error")
    ylabel("Bit Error Rate")
    title("Incremental Channel Coding BER")

    %plotting the Throughput values
    figure()
    plot(probz,Throughput)
    xlabel("Probability of Error")
    ylabel("Throughput")
    title("Incremental Channel Coding Throughput")
