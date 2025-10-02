function Y = NR_LDPC_client(server_ip, server_port, data, length_data)
  uByte = udpport("byte", "Timeout", 10); % Increase timeout to 10 seconds

  write(uByte, data, "double", server_ip, server_port);

  try
    Y = read(uByte, length_data -1, "uint8");
  catch ME
    warning("Error during data transmission: %s", MEx.message);
    Y = [];
  end

  flush(uByte)
