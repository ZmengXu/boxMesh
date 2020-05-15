/* Creates the same files that would be created by OpenFOAM's blockMesh if the
 * mesh is rectangular cuboidal domain divided into smaller identical
 * rectangular cuboids.
 * Using this code, I was able to generate a 500x500x500 mesh that I couldn't
 * generate using blockMesh on my computer with 64 GB RAM. It takes 70% more
 * time than the OpenFOAM's blockMesh, but requires 40% less memory.
 *
 * See comments on void boxMesh(int Nx, int Ny, int Nz, double Lx, double Ly, double Lz)
 * below for an example of the equivalent blockMeshDict dictionary
 *
 * This code requires that the directory constant/polyMesh exits.
 *
 * Author: Thomas Oliveira
 * Contact: Gmail - thomas.oliveira */

#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <vector>

#include <fstream>
#include <iostream>
#include <sys/stat.h>

using namespace std;

/* defines how many iterations are done before ostringstream is streamed to
 * ofstream and cleared. Clearing ostringstream saves memory */
#define CLEAR_OSTRINGSTREAM 1000000

typedef vector<int> Points;

vector<Points> pointsOfInternalFacesOfCell(int cellId, int Nx, int Ny, int Nz);
vector<int> ownersOfInternalFacesOfCell(int cellId, int Nx, int Ny, int Nz);
vector<int> neighboursOfInternalFacesOfCell(int cellId, int Nx, int Ny, int Nz);
Points pointsOfExternalFace(int externalFaceId, int Nx, int Ny, int Nz);
int ownerOfExternalFace(int externalFaceId, int Nx, int Ny, int Nz);

void writeBoundary(int nFaces, int nExternalFaces)
{
    ofstream currentFile("constant/polyMesh/boundary");

    currentFile << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    currentFile << "| =========                 |                                                 |\n";
    currentFile << "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n";
    currentFile << "|  \\\\    /   O peration     | Version:  2.3.0                                 |\n";
    currentFile << "|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n";
    currentFile << "|    \\\\/     M anipulation  |                                                 |\n";
    currentFile << "\\*---------------------------------------------------------------------------*/\n";
    currentFile << "FoamFile\n";
    currentFile << "{\n";
    currentFile << "    version     2.0;\n";
    currentFile << "    format      ascii;\n";
    currentFile << "    class       polyBoundaryMesh;\n";
    currentFile << "    location    \"constant/polyMesh\";\n";
    currentFile << "    object      boundary;\n";
    currentFile << "}\n";
    currentFile << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";
    currentFile << "\n";
    currentFile << "1\n";
    currentFile << "(\n";
    currentFile << "    defaultFaces\n";
    currentFile << "    {\n";
    currentFile << "        type            empty;\n";
    currentFile << "        inGroups        1(empty);\n";

    ostringstream convert;
    convert << "        nFaces          " << nExternalFaces        << ";\n";
    convert << "        startFace       " << nFaces-nExternalFaces << ";\n";
    currentFile << convert.str();

    currentFile << "    }\n";
    currentFile << ")\n";
    currentFile << "\n";
    currentFile << "// ************************************************************************* //\n";
    currentFile.close();
}

void writeFaces(int nCells, int nFaces, int nExternalFaces, int Nx, int Ny, int Nz)
{
    ofstream currentFile("constant/polyMesh/faces");

    currentFile << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    currentFile << "| =========                 |                                                 |\n";
    currentFile << "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n";
    currentFile << "|  \\\\    /   O peration     | Version:  2.3.0                                 |\n";
    currentFile << "|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n";
    currentFile << "|    \\\\/     M anipulation  |                                                 |\n";
    currentFile << "\\*---------------------------------------------------------------------------*/\n";
    currentFile << "FoamFile\n";
    currentFile << "{\n";
    currentFile << "    version     2.0;\n";
    currentFile << "    format      ascii;\n";
    currentFile << "    class       faceList;\n";
    currentFile << "    location    \"constant/polyMesh\";\n";
    currentFile << "    object      faces;\n";
    currentFile << "}\n";
    currentFile << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";
    currentFile << "\n";
    currentFile << "\n";

    ostringstream convert;
    convert << nFaces << "\n(\n";

    vector<int> points;
    vector<Points> faces;
    int size;

    int countFaces = 0;
    for (int i=0; i<nCells; i++)
    {
        faces = pointsOfInternalFacesOfCell(i, Nx, Ny, Nz);

        size = faces.size();
        for (int j=0; j<size; j++)
        {
            convert << "4(" << faces[j][0] << " " << faces[j][1] << " " << faces[j][2] << " " << faces[j][3] << ")\n";
            countFaces++;
        }
        if (countFaces>CLEAR_OSTRINGSTREAM)
        {
            countFaces = 0;
            currentFile << convert.str();
            convert.clear();
            convert.str("");
        }
    }

    countFaces = 0;
    for (int i=0; i<nExternalFaces; i++)
    {
        points = pointsOfExternalFace(i, Nx, Ny, Nz);
        convert << "4(" <<  points[0] << " " << points[1] << " " << points[2] << " " << points[3] << ")\n";
        countFaces++;
        if (countFaces>CLEAR_OSTRINGSTREAM)
        {
            countFaces = 0;
            currentFile << convert.str();
            convert.clear();
            convert.str("");
        }
    }
    currentFile << convert.str();

    currentFile << ")\n";
    currentFile << "\n";
    currentFile << "\n";
    currentFile << "// ************************************************************************* //\n";
    currentFile.close();
}

void writeNeighbour(int nPoints, int nCells, int nFaces, int nExternalFaces, int Nx, int Ny, int Nz)
{
    ofstream currentFile("constant/polyMesh/neighbour");

    currentFile << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    currentFile << "| =========                 |                                                 |\n";
    currentFile << "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n";
    currentFile << "|  \\\\    /   O peration     | Version:  2.3.0                                 |\n";
    currentFile << "|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n";
    currentFile << "|    \\\\/     M anipulation  |                                                 |\n";
    currentFile << "\\*---------------------------------------------------------------------------*/\n";
    currentFile << "FoamFile\n";
    currentFile << "{\n";
    currentFile << "    version     2.0;\n";
    currentFile << "    format      ascii;\n";
    currentFile << "    class       labelList;\n";

    ostringstream convert;
    convert << "    note        \"nPoints: " << nPoints << " nCells: " << nCells << " nFaces: " << nFaces << " nInternalFaces: " << nFaces-nExternalFaces << "\";\n";
    currentFile << convert.str();

    currentFile << "    location    \"constant/polyMesh\";\n";
    currentFile << "    object      neighbour;\n";
    currentFile << "}\n";
    currentFile << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";
    currentFile << "\n";
    currentFile << "\n";

    convert.clear();
    convert.str("");
    convert << nFaces-nExternalFaces << "\n(\n";

    vector<int> neighbours;
    int size;

    int countFaces = 0;
    for (int i=0; i<nCells; i++)
    {
        neighbours = neighboursOfInternalFacesOfCell(i, Nx, Ny, Nz);

        size = neighbours.size();
        for (int j=0; j<size; j++)
        {
            convert << neighbours[j] << "\n";
            countFaces++;
        }
        if (countFaces>CLEAR_OSTRINGSTREAM)
        {
            countFaces = 0;
            currentFile << convert.str();
            convert.clear();
            convert.str("");
        }
    }

    currentFile << convert.str();

    currentFile << ")\n";
    currentFile << "\n";
    currentFile << "\n";
    currentFile << "// ************************************************************************* //\n";
    currentFile.close();
}

void writeOwner(int nPoints, int nCells, int nFaces, int nExternalFaces, int Nx, int Ny, int Nz)
{
    ofstream currentFile("constant/polyMesh/owner");

    currentFile << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    currentFile << "| =========                 |                                                 |\n";
    currentFile << "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n";
    currentFile << "|  \\\\    /   O peration     | Version:  2.3.0                                 |\n";
    currentFile << "|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n";
    currentFile << "|    \\\\/     M anipulation  |                                                 |\n";
    currentFile << "\\*---------------------------------------------------------------------------*/\n";
    currentFile << "FoamFile\n";
    currentFile << "{\n";
    currentFile << "    version     2.0;\n";
    currentFile << "    format      ascii;\n";
    currentFile << "    class       labelList;\n";

    ostringstream convert;
    convert << "    note        \"nPoints: " << nPoints << " nCells: " << nCells << " nFaces: " << nFaces << " nInternalFaces: " << nFaces-nExternalFaces << "\";\n";
    currentFile << convert.str();

    currentFile << "    location    \"constant/polyMesh\";\n";
    currentFile << "    object      owner;\n";
    currentFile << "}\n";
    currentFile << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";
    currentFile << "\n";
    currentFile << "\n";

    convert.clear();
    convert.str("");
    convert << nFaces << "\n(\n";

    vector<int> owners;
    int size;

    int countFaces = 0;
    for (int i=0; i<nCells; i++)
    {
        owners = ownersOfInternalFacesOfCell(i, Nx, Ny, Nz);

        size = owners.size();
        for (int j=0; j<size; j++)
        {
            convert << owners[j] << "\n";
            countFaces++;
        }
        if (countFaces>CLEAR_OSTRINGSTREAM)
        {
            countFaces = 0;
            currentFile << convert.str();
            convert.clear();
            convert.str("");
        }
    }

    countFaces = 0;
    for (int i=0; i<nExternalFaces; i++)
    {
        convert << ownerOfExternalFace(i, Nx, Ny, Nz)  << "\n";
        countFaces++;
        if (countFaces>CLEAR_OSTRINGSTREAM)
        {
            countFaces = 0;
            currentFile << convert.str();
            convert.clear();
            convert.str("");
        }
    }
    currentFile << convert.str();

    currentFile << ")\n";
    currentFile << "\n";
    currentFile << "\n";
    currentFile << "// ************************************************************************* //\n";
    currentFile.close();
}

void writePoints(int nPoints, int Nx, int Ny, int Nz, double Lx, double Ly, double Lz)
{
    ofstream currentFile("constant/polyMesh/points");

    currentFile << "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    currentFile << "| =========                 |                                                 |\n";
    currentFile << "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n";
    currentFile << "|  \\\\    /   O peration     | Version:  2.3.0                                 |\n";
    currentFile << "|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n";
    currentFile << "|    \\\\/     M anipulation  |                                                 |\n";
    currentFile << "\\*---------------------------------------------------------------------------*/\n";
    currentFile << "FoamFile\n";
    currentFile << "{\n";
    currentFile << "    version     2.0;\n";
    currentFile << "    format      ascii;\n";
    currentFile << "    class       vectorField;\n";
    currentFile << "    location    \"constant/polyMesh\";\n";
    currentFile << "    object      points;\n";
    currentFile << "}\n";
    currentFile << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";
    currentFile << "\n";
    currentFile << "\n";

    ostringstream convert;
    convert.clear();
    convert.str("");
    convert << nPoints << "\n(\n";
    currentFile << convert.str();

    int xi, zi, yi;
    for (int i=0; i<nPoints; i++)
    {
        // Calculate point grid coordinates
        xi = i % (Nx+1);
        zi = int(i/((Nx+1)*(Ny+1)));
        yi = ((i-(Nx+1)*(Ny+1)*zi)-xi)/(Nx+1);

        convert.clear();
        convert.str("");
        //  Write point spatial coordinates
        convert << "(" << setprecision(10) << xi*Lx << " " << yi*Ly << " " << zi*Lz << ")\n";
        currentFile << convert.str();
    }

    currentFile << ")\n";
    currentFile << "\n";
    currentFile << "\n";
    currentFile << "// ************************************************************************* //\n";
    currentFile.close();
}

vector<Points> pointsOfInternalFacesOfCell(int cellId, int Nx, int Ny, int Nz)
{
    /* Returns a vector of Points of cell internal faces with outward normals
     * (1,0,0), (0,1,0), (0,0,1), in this order. */

    vector<Points> output;
    Points points;
    int point_1, point_2, point_3, point_4;

    // grid coordinates of the cell
    int xCell = cellId % Nx;
    int zCell = int(cellId/(Nx*Ny));
    int yCell = ((cellId-Nx*Ny*zCell)-xCell)/Nx;

    if (xCell < Nx-1) // Check if face with outward normal (1,0,0) is internal
    {
        point_1 = zCell*(Nx+1)*(Ny+1) + yCell*(Nx+1) + xCell + 1;
        point_2 = point_1 + (Nx+1);
        point_3 = point_2 + (Nx+1)*(Ny+1);
        point_4 = point_1 + (Nx+1)*(Ny+1);

        points.clear();
        points.push_back(point_1);
        points.push_back(point_2);
        points.push_back(point_3);
        points.push_back(point_4);

        output.push_back(points);
    }

    if (yCell < Ny-1) // Check if face with outward normal (0,1,0) is internal
    {
        point_1 = zCell*(Nx+1)*(Ny+1) + (yCell+1)*(Nx+1) + xCell;
        point_2 = point_1 + (Nx+1)*(Ny+1);
        point_3 = point_2 + 1;
        point_4 = point_1 + 1;

        points.clear();
        points.push_back(point_1);
        points.push_back(point_2);
        points.push_back(point_3);
        points.push_back(point_4);

        output.push_back(points);
    }

    if (zCell < Nz-1) // Check if face with outward normal (0,0,1) is internal
    {
        point_1 = (zCell+1)*(Nx+1)*(Ny+1) + yCell*(Nx+1) + xCell;
        point_2 = point_1 + 1;
        point_3 = point_2 + (Nx+1);
        point_4 = point_1 + (Nx+1);

        points.clear();
        points.push_back(point_1);
        points.push_back(point_2);
        points.push_back(point_3);
        points.push_back(point_4);

        output.push_back(points);
    }

    return output;
}

vector<int> ownersOfInternalFacesOfCell(int cellId, int Nx, int Ny, int Nz)
{
    /* Returns a vector of owners of cell internal faces with outward normals
     * (1,0,0), (0,1,0), (0,0,1), in this order. */

    vector<int> output;
    int owner;

    // grid coordinates of the cell
    int xCell = cellId % Nx;
    int zCell = int(cellId/(Nx*Ny));
    int yCell = ((cellId-Nx*Ny*zCell)-xCell)/Nx;

    if (xCell < Nx-1) // Check if face with outward normal (1,0,0) is internal
    {
        owner = cellId;
        output.push_back(owner);
    }

    if (yCell < Ny-1) // Check if face with outward normal (0,1,0) is internal
    {
        owner = cellId;
        output.push_back(owner);
    }

    if (zCell < Nz-1) // Check if face with outward normal (0,0,1) is internal
    {
        owner = cellId;
        output.push_back(owner);
    }

    return output;
}

vector<int> neighboursOfInternalFacesOfCell(int cellId, int Nx, int Ny, int Nz)
{
    /* Returns a vector of neighbours of cell internal faces with outward normals
     * (1,0,0), (0,1,0), (0,0,1), in this order. */

    vector<int> output;
    int neighbour;

    // grid coordinates of the cell
    int xCell = cellId % Nx;
    int zCell = int(cellId/(Nx*Ny));
    int yCell = ((cellId-Nx*Ny*zCell)-xCell)/Nx;

    if (xCell < Nx-1) // Check if face with outward normal (1,0,0) is internal
    {
        neighbour = cellId + 1;
        output.push_back(neighbour);
    }

    if (yCell < Ny-1) // Check if face with outward normal (0,1,0) is internal
    {
        neighbour = cellId + Nx;
        output.push_back(neighbour);
    }

    if (zCell < Nz-1) // Check if face with outward normal (0,0,1) is internal
    {
        neighbour = cellId + Nx*Ny;
        output.push_back(neighbour);
    }

    return output;
}

Points pointsOfExternalFace(int externalFaceId, int Nx, int Ny, int Nz)
{
    /* Returns Points of face nInternalFaces+externalFaceId (which is an external face) */

    /* The way to calculate the points depends on which face of the domain it belongs to:
     * The first Ny*Nz faces have outward normal (-1, 0, 0) (referred to as x-)
     * The next  Ny*Nz faces have outward normal ( 1, 0, 0) (referred to as x+)
     * The next  Nx*Nz faces have outward normal ( 0,-1, 0) (referred to as y-)
     * The next  Nx*Nz faces have outward normal ( 0, 1, 0) (referred to as y+)
     * The next  Nx*Ny faces have outward normal ( 0, 0,-1) (referred to as z-)
     * The next  Nx*Ny faces have outward normal ( 0, 0, 1) (referred to as z+)
     * The following tests check in which face of the domain the face is. */
    string normal;
    if   (externalFaceId < Ny*Nz)
        normal = "x-";
    else if (externalFaceId < 2*Ny*Nz)
        normal = "x+";
    else if (externalFaceId < 2*Ny*Nz + Nx*Nz)
        normal = "y-";
    else if (externalFaceId < 2*Ny*Nz + 2*Nx*Nz)
        normal = "y+";
    else if (externalFaceId < 2*Ny*Nz + 2*Nx*Nz + Nx*Ny)
        normal = "z-";
    else if (externalFaceId < 2*Ny*Nz + 2*Nx*Nz + 2*Nx*Ny)
        normal = "z+";
    else
        throw invalid_argument("externalFaceId is too large");

    /* Now that the normal to the face is known, calculate the points in the
     * respective way. */
    int temp;
    int point_1, point_2, point_3, point_4;
    if (normal == "x-")
    {
        temp  = externalFaceId;     // index on "x-" face
        int zCell = int(temp/Ny);   // z-grid-coordinate of face
        int yCell = temp%Ny;        // y-grid-coordinate of face
        point_1 = zCell*(Nx+1)*(Ny+1) + yCell*(Nx+1);
        point_2 = point_1 + (Nx+1)*(Ny+1);
        point_3 = point_2 + (Nx+1);
        point_4 = point_1 + (Nx+1);
    }
    else if (normal == "x+")
    {
        /* To get the formula for point_1, start with the formula for point_1 of face x-,
         * and replace externalFaceId by (externalFaceId-Ny*Nz), and then add Nx. */
        temp  = externalFaceId - (Ny*Nz);
        int zCell = int(temp/Ny);
        int yCell = temp%Ny;
        point_1 = zCell*(Nx+1)*(Ny+1) + yCell*(Nx+1) + Nx;
        point_2 = point_1 + (Nx+1);
        point_3 = point_2 + (Nx+1)*(Ny+1);
        point_4 = point_1 + (Nx+1)*(Ny+1);
    }
    else if (normal == "y-")
    {
        temp  = externalFaceId - (2*Ny*Nz);     // index on "y-" face
        int xCell = int(temp/Nz);               // x-grid-coordinate of face
        int zCell = temp%Nz;                    // z-grid-coordinate of face
        point_1 = zCell*(Nx+1)*(Ny+1) + xCell;
        point_2 = point_1 + 1;
        point_3 = point_2 + (Nx+1)*(Ny+1);
        point_4 = point_3 -1;
    }
    else if (normal == "y+")
    {
        /* To get the formula for point_1, start with the formula for point_1 of face y-,
         * and replace externalFaceId by (externalFaceId-Nx*Nz), and then add (Nx+1)*Ny. */
        temp  = externalFaceId - (2*Ny*Nz + Nx*Nz);
        int xCell = int(temp/Nz);
        int zCell = temp%Nz;
        point_1 = zCell*(Nx+1)*(Ny+1) + xCell + (Nx+1)*Ny;
        point_2 = point_1 + (Nx+1)*(Ny+1);
        point_3 = point_2 + 1;
        point_4 = point_1 + 1;
    }
    else if (normal == "z-")
    {
        temp  = externalFaceId - (2*Ny*Nz + 2*Nx*Nz);   // index on "z-" face
        int xCell = int(temp/Ny);                       // x-grid-coordinate of face
        int yCell = temp%Ny;                            // y-grid-coordinate of face
        point_1 = yCell*(Nx+1) + xCell;
        point_2 = point_1 + (Nx+1);
        point_3 = point_2 + 1;
        point_4 = point_1 + 1;
    }
    else if (normal == "z+")
    {
        /* To get the formula for point_1, start with the formula for point_1 of face z-,
         * and replace externalFaceId by (externalFaceId-Nx*Ny), and then add (Nx+1)*(Ny+1)*Nz. */
        temp  = externalFaceId - (2*Ny*Nz + 2*Nx*Nz + Nx*Ny);
        int xCell = int(temp/Ny);
        int yCell = temp%Ny;
        point_1 = yCell*(Nx+1) + xCell + (Nx+1)*(Ny+1)*Nz;
        point_2 = point_1 + 1;
        point_3 = point_2 + (Nx+1);
        point_4 = point_3 -1;
    }
    else
        throw logic_error("Unexpected value for normal: " + normal);

    Points points;
    points.push_back(point_1);
    points.push_back(point_2);
    points.push_back(point_3);
    points.push_back(point_4);

    return points;
}

int ownerOfExternalFace(int externalFaceId, int Nx, int Ny, int Nz)
{
    /* Returns owner of face nInternalFaces+externalFaceId (which is an external face) */

    /* The way to calculate the owner depends on which face of the domain it belongs to:
     * The first Ny*Nz faces have outward normal (-1, 0, 0) (referred to as x-)
     * The next  Ny*Nz faces have outward normal ( 1, 0, 0) (referred to as x+)
     * The next  Nx*Nz faces have outward normal ( 0,-1, 0) (referred to as y-)
     * The next  Nx*Nz faces have outward normal ( 0, 1, 0) (referred to as y+)
     * The next  Nx*Ny faces have outward normal ( 0, 0,-1) (referred to as z-)
     * The next  Nx*Ny faces have outward normal ( 0, 0, 1) (referred to as z+)
     * The following tests check in which face of the domain the face is. */
    string normal;
    if   (externalFaceId < Ny*Nz)
        normal = "x-";
    else if (externalFaceId < 2*Ny*Nz)
        normal = "x+";
    else if (externalFaceId < 2*Ny*Nz + Nx*Nz)
        normal = "y-";
    else if (externalFaceId < 2*Ny*Nz + 2*Nx*Nz)
        normal = "y+";
    else if (externalFaceId < 2*Ny*Nz + 2*Nx*Nz + Nx*Ny)
        normal = "z-";
    else if (externalFaceId < 2*Ny*Nz + 2*Nx*Nz + 2*Nx*Ny)
        normal = "z+";
    else
        throw invalid_argument("externalFaceId is too large");

    /* Now that the normal to the face is known, calculate the owner in the
     * respective way. */
    int temp;
    int owner;
    if (normal == "x-")
    {
        temp  = externalFaceId;     // index on "x-" face
        int zCell = int(temp/Ny);   // z-grid-coordinate of face
        int yCell = temp%Ny;        // y-grid-coordinate of face
        owner = yCell*Nx + zCell*(Nx*Ny);
    }
    else if (normal == "x+")
    {
        temp  = externalFaceId - (Ny*Nz);
        int zCell = int(temp/Ny);
        int yCell = temp%Ny;
        owner = yCell*Nx + zCell*(Nx*Ny) + (Nx-1);
    }
    else if (normal == "y-")
    {
        temp  = externalFaceId - (2*Ny*Nz);     // index on "y-" face
        int xCell = int(temp/Nz);               // x-grid-coordinate of face
        int zCell = temp%Nz;                    // z-grid-coordinate of face
        owner = xCell + zCell*(Nx*Ny);
    }
    else if (normal == "y+")
    {
        temp  = externalFaceId - (2*Ny*Nz + Nx*Nz);
        int xCell = int(temp/Nz);
        int zCell = temp%Nz;
        owner = xCell + zCell*(Nx*Ny) + Nx*(Ny-1);
    }
    else if (normal == "z-")
    {
        temp  = externalFaceId - (2*Ny*Nz + 2*Nx*Nz);   // index on "z-" face
        int xCell = int(temp/Ny);                       // x-grid-coordinate of face
        int yCell = temp%Ny;                            // y-grid-coordinate of face
        owner = xCell + yCell*Nx;
    }
    else if (normal == "z+")
    {
        temp  = externalFaceId - (2*Ny*Nz + 2*Nx*Nz + Nx*Ny);
        int xCell = int(temp/Ny);
        int yCell = temp%Ny;
        owner = xCell + yCell*Nx + Nx*Ny*(Nz-1);
    }
    else
        throw logic_error("Unexpected value for normal: " + normal);

    return owner;
}

void boxMesh(int Nx, int Ny, int Nz, double Lx, double Ly, double Lz)
{
    /* Create files
        constant/polyMesh/boundary
        constant/polyMesh/faces
        constant/polyMesh/neighbour
        constant/polyMesh/owner
        constant/polyMesh/points
    that would be created by calling OpenFOAM's blockMesh using the following
    blockMeshDict dictionary:

        vertices
        (
            (0     0     0    )
            (Nx*Lx 0     0    )
            (Nx*Lx Ny*Ly 0    )
            (0     Ny*Ly 0    )
            (0     0     Nz*Lz)
            (Nx*Lx 0     Nz*Lz)
            (Nx*Lx Ny*Ly Nz*Lz)
            (0     Ny*Ly Nz*Lz)
        );

        blocks
        (
            hex (0 1 2 3 4 5 6 7) (Nx Ny Nz) simpleGrading (1 1 1)
        );

        edges
        (
        );

        boundary
        (
        );

        mergePatchPairs
        (
        );
    */

    // Calculate number of points, faces and cells
    int nPoints        = (Nx+1)*(Ny+1)*(Nz+1);
    int nCells         = Nx*Ny*Nz;
    int nFaces         = Nx*Ny*(Nz+1)+Nx*Nz*(Ny+1)+Ny*Nz*(Nx+1);
    int nExternalFaces = 2*(Nx*Ny+Nx*Nz+Ny*Nz);

    // Write files
    writeBoundary(nFaces, nExternalFaces);
    writeFaces(nCells, nFaces, nExternalFaces, Nx, Ny, Nz);
    writeNeighbour(nPoints, nCells, nFaces, nExternalFaces, Nx, Ny, Nz);
    writeOwner(nPoints, nCells, nFaces, nExternalFaces, Nx, Ny, Nz);
    writePoints(nPoints, Nx, Ny, Nz, Lx, Ly, Lz);
}

int main(int argc, char *argv[])
{
    int Nx, Ny, Nz;
    double Lx, Ly, Lz;
    if ( argc != 7 ) // argc should be 2 for correct execution
        // We print argv[0] assuming it is the program name
        cout << "usage: " << argv[0] << " Nx Ny Nz Lx Ly Lz\n";
    else
    {
        Nx = atoi(argv[1]);
        Ny = atoi(argv[2]);
        Nz = atoi(argv[3]);
        Lx = atof(argv[4]);
        Ly = atof(argv[5]);
        Lz = atof(argv[6]);

        struct stat info;
        if( stat("constant/polyMesh", &info)==0 )   // Check if directory exists.
        boxMesh(Nx,Ny,Nz,Lx,Ly,Lz);
        else
            cout << "Please create directory constant/polyMesh first\n";
    }
    return 0;
}

