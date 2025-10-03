function Y = NR_LDPC_client(server_ip, server_port, data, length_data)
  uByte = udpport("byte", "Timeout", 10); % Increase timeout to 10 seconds
  
  write(uByte, data, "double", server_ip, server_port);

  try
    Y = read(uByte, length_data, "uint8");
  catch MException
    warning("Error during data transmission!");
    Y = [];
  end

  flush(uByte)
