// : SocketHandler.h
#pragma once
#include <string>

bool InitializeSocket(int port);
void SendData(const std::string& data, const std::string& client_ip, int client_port);
void CleanupSocket();
