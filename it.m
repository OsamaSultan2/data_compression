clc; 
clear all; 
close all;

%% Ask user to select an image
 [filename, pathname] = uigetfile({'*.jpg;*.jpeg;*.png;*.bmp;*.tif;*.webp', 'Image Files'}, 'Select an Image');
 if isequal(filename,0) || isequal(pathname,0)
     disp('No image selected. Exiting.');
     return;
 else
     selected_file = fullfile(pathname, filename);
     disp(['User selected ', selected_file])
 end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load and process the image
  img_path = selected_file;
  original_image = imread(img_path);

% Convert to grayscale if needed
  if size(original_image, 3) == 3
      original_image = rgb2gray(original_image);
  end
    
% Display original image
  figure('Name', 'Original Image', 'NumberTitle', 'off');
  imshow(original_image);
  title('Original Grayscale Image');

% Get image dimensions
  [rows, cols] = size(original_image);
  fprintf('Image dimensions: %d x %d\n', rows, cols);
  fprintf('Total pixels: %d\n', rows * cols);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load and process the image
   data = reshape(original_image.', 1, []);    % Stack rows into a single row

% Calculate probabilities of each symbol 
  symbols = unique(data); % Remove repeated symbols 
  probabilities = zeros(size(symbols)); % Initialize probability array 
  
  for i = 1:length(symbols) 
      probabilities(i) = sum(data == symbols(i)) / length(data); 
  end 

% Calculate information for each symbol 
  info = -log2(probabilities); 
  info_to_ent = probabilities .* info; 

% Calculate entropy 
  entropy = sum(info_to_ent(isfinite(info_to_ent))); 


% Display information and entropy
  fprintf('%-12s | %-12s | %-16s\n','Gray level(Symbols)','Probability','Information gained');
  fprintf('---------------------------------------------------------------\n');

  for i = 1:length(symbols)
      fprintf('%-19d | %-12.5f | %-16.4f\n', symbols(i), probabilities(i), info(i));
  end

  fprintf('---------------------------------------------------------------\n');
  fprintf('Entropy: %.6f bits/symbol\n', entropy);

%% Huffman Coding 
codeword_length = length(probabilities); 
codeword = cell(1, codeword_length); % Empty cell array for Huffman codewords 
x = zeros(codeword_length, codeword_length); % Matrix for Huffman algorithm 
% Temporary probability array to keep original values intact 
temp_p = probabilities; 

% Huffman algorithm 
for a = 1:codeword_length-1 
[~, idx] = sort(temp_p); % Sort probabilities 
% Combine the smallest two probabilities 
temp_p(idx(2)) = temp_p(idx(1)) + temp_p(idx(2));  
x(idx(1), a) = 11; % Assign value to the first smallest p 
x(idx(2), a) = 10; % Assign value to the second smallest p 
temp_p(idx(1)) = nan; % Mark the combined node 
end 
% Generate Huffman codewords based on the matrix 
for a = codeword_length-1:-1:1 
code_1 = find(x(:, a) == 11); 
code_0 = find(x(:, a) == 10); 
codeword{code_1} = [codeword{code_0}, '1']; 
codeword{code_0} = [codeword{code_0}, '0']; 
end 
% Display Huffman dictionary with ASCII code 
disp('---------------------------------'); 
disp('Huffman Coding: Symbol | Codeword'); 
for i = 1:length(symbols) 
fprintf('%c | %d | %s\n', symbols(i), double(symbols(i)), codeword{i}); 
end 
% Encode the data using Huffman codes 
encoded = strings(1, length(data)); 
for i = 1:length(data) 
for j = 1:length(symbols) 
if data(i) == symbols(j) 
encoded(i) = codeword{j}; 
end 
end 
end 
% Transmitted Binary Code (Encoded Message - Huffman) 
str = strjoin(encoded); 
str = strrep(str, ' ', ''); 
disp('---------------------------------'); 
disp('Transmitted Binary Code (Encoded Message - Huffman):'); 
disp(str); 
% Decode the message using Huffman codes 
decoded = []; % Initialize decoded data array 
stream = []; % Initialize temporary stream 
for i = 1:length(encoded) 
stream = [stream, encoded(i)]; 
idx = find(strcmp(stream, codeword)); 
if ~isempty(idx) 
decoded = [decoded, symbols(idx)]; 
stream = []; 
end 
end 
disp('---------------------------------'); 
disp('Decoded Data (Decoded Message - Huffman):'); 
%disp(decoded); 