// : SocketHandler.cpp
#include "SocketHandler.h"
#include <WinSock2.h>
#include <iostream>
#include <ws2tcpip.h>

#pragma comment(lib, "ws2_32.lib")

SOCKET server_socket;

bool InitializeSocket(int port) {
    WSADATA wsaData;
    int iResult = WSAStartup(MAKEWORD(2, 2), &wsaData);
    if (iResult != NO_ERROR) {
        std::cerr << "WSAStartup failed with error: " << iResult << std::endl;
        return false;
    }

    server_socket = socket(AF_INET, SOCK_DGRAM, IPPROTO_UDP);
    if (server_socket == INVALID_SOCKET) {
        std::cerr << "Socket creation failed with error: " << WSAGetLastError() << std::endl;
        WSACleanup();
        return false;
    }

    sockaddr_in server_address;
    server_address.sin_family = AF_INET;
    server_address.sin_port = htons(port);
    server_address.sin_addr.S_un.S_addr = INADDR_ANY;

    if (bind(server_socket, (SOCKADDR*)&server_address, sizeof(server_address)) == SOCKET_ERROR) {
        std::cerr << "Bind failed with error: " << WSAGetLastError() << std::endl;
        closesocket(server_socket);
        WSACleanup();
        return false;
    }

    return true;
}

void SendData(const std::string& data, const std::string& client_ip, int client_port) {
    sockaddr_in client_address;
    client_address.sin_family = AF_INET;
    client_address.sin_port = htons(client_port);
    inet_pton(AF_INET, client_ip.c_str(), &client_address.sin_addr);

    int sendResult = sendto(server_socket, data.c_str(), data.length(), 0, (SOCKADDR*)&client_address, sizeof(client_address));
    if (sendResult == SOCKET_ERROR) {
        std::cerr << "Send failed with error: " << WSAGetLastError() << std::endl;
    }
}

void CleanupSocket() {
    closesocket(server_socket);
    WSACleanup();
}
