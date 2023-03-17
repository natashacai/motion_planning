/*
    MIT License

    Copyright (c) 2021 Zhepei Wang (wangzhepei@live.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP

#include <Eigen/Eigen>

#include <iostream>
#include <cmath>
#include <cfloat>
#include <vector>

template <int D>
class Piece
{
public:
    typedef Eigen::Matrix<double, 3, D + 1> CoefficientMat;
    typedef Eigen::Matrix<double, 3, D> VelCoefficientMat;
    typedef Eigen::Matrix<double, 3, D - 1> AccCoefficientMat;

private:
    double duration;
    CoefficientMat coeffMat;

public:
    Piece() = default;

    Piece(double dur, const CoefficientMat &cMat)
        : duration(dur), coeffMat(cMat) {}

    inline int getDim() const
    {
        return 3;
    }

    inline int getDegree() const
    {
        return D;
    }

    inline double getDuration() const
    {
        return duration;
    }

    inline const CoefficientMat &getCoeffMat() const
    {
        return coeffMat;
    }

    inline Eigen::Vector3d getPos(const double &t) const
    {//get position
        Eigen::Vector3d pos(0.0, 0.0, 0.0);
        double tn = 1.0;
        for (int i = D; i >= 0; i--)
        {
            pos += tn * coeffMat.col(i);
            tn *= t;
        }
        return pos;
    }

    inline Eigen::Vector3d getVel(const double &t) const
    { // get velocity
        Eigen::Vector3d vel(0.0, 0.0, 0.0);
        double tn = 1.0;
        int n = 1;
        for (int i = D - 1; i >= 0; i--)
        {
            vel += n * tn * coeffMat.col(i);
            tn *= t;
            n++;
        }
        return vel;
    }

    inline Eigen::Vector3d getAcc(const double &t) const
    {//get acceleration
        Eigen::Vector3d acc(0.0, 0.0, 0.0);
        double tn = 1.0;
        int m = 1;
        int n = 2;
        for (int i = D - 2; i >= 0; i--)
        {
            acc += m * n * tn * coeffMat.col(i);
            tn *= t;
            m++;
            n++;
        }
        return acc;
    }

    inline Eigen::Vector3d getJer(const double &t) const
    { //get jerk
        Eigen::Vector3d jer(0.0, 0.0, 0.0);
        double tn = 1.0;
        int l = 1;
        int m = 2;
        int n = 3;
        for (int i = D - 3; i >= 0; i--)
        {
            jer += l * m * n * tn * coeffMat.col(i);
            tn *= t;
            l++;
            m++;
            n++;
        }
        return jer;
    }

    inline CoefficientMat normalizePosCoeffMat() const
    { // normalized position coefficients
        CoefficientMat nPosCoeffsMat;
        double t = 1.0;
        for (int i = D; i >= 0; i--)
        {
            nPosCoeffsMat.col(i) = coeffMat.col(i) * t;
            t *= duration;
        }
        return nPosCoeffsMat;
    }

    inline VelCoefficientMat normalizeVelCoeffMat() const
    { // normalized velocity coefficients
        VelCoefficientMat nVelCoeffMat;
        int n = 1;
        double t = duration;
        for (int i = D - 1; i >= 0; i--)
        {
            nVelCoeffMat.col(i) = n * coeffMat.col(i) * t;
            t *= duration;
            n++;
        }
        return nVelCoeffMat;
    }

    inline AccCoefficientMat normalizeAccCoeffMat() const
    { // normalized accelaration coefficients
        AccCoefficientMat nAccCoeffMat;
        int n = 2;
        int m = 1;
        double t = duration * duration;
        for (int i = D - 2; i >= 0; i--)
        {
            nAccCoeffMat.col(i) = n * m * coeffMat.col(i) * t;
            n++;
            m++;
            t *= duration;
        }
        return nAccCoeffMat;
    }
};

template <int D>
class Trajectory
{
private:
    typedef std::vector<Piece<D>> Pieces;
    Pieces pieces;

public:
    Trajectory() = default;

    Trajectory(const std::vector<double> &durs,
               const std::vector<typename Piece<D>::CoefficientMat> &cMats)
    {
        int N = std::min(durs.size(), cMats.size());
        pieces.reserve(N);
        for (int i = 0; i < N; i++)
        {
            pieces.emplace_back(durs[i], cMats[i]);
        }
    }

    inline int getPieceNum() const
    {
        return pieces.size();
    }

    inline Eigen::VectorXd getDurations() const
    {
        int N = getPieceNum();
        Eigen::VectorXd durations(N);
        for (int i = 0; i < N; i++)
        {
            durations(i) = pieces[i].getDuration();
        }
        return durations;
    }

    inline double getTotalDuration() const
    {
        int N = getPieceNum();
        double totalDuration = 0.0;
        for (int i = 0; i < N; i++)
        {
            totalDuration += pieces[i].getDuration();
        }
        return totalDuration;
    }

    inline Eigen::Matrix3Xd getPositions() const
    {
        int N = getPieceNum();
        Eigen::Matrix3Xd positions(3, N + 1);
        for (int i = 0; i < N; i++)
        {
            positions.col(i) = pieces[i].getCoeffMat().col(D);
        }
        positions.col(N) = pieces[N - 1].getPos(pieces[N - 1].getDuration());
        return positions;
    }

    inline const Piece<D> &operator[](int i) const
    {
        return pieces[i];
    }

    inline Piece<D> &operator[](int i)
    {
        return pieces[i];
    }

    inline void clear(void)
    {
        pieces.clear();
        return;
    }

    inline typename Pieces::const_iterator begin() const
    {
        return pieces.begin();
    }

    inline typename Pieces::const_iterator end() const
    {
        return pieces.end();
    }

    inline typename Pieces::iterator begin()
    {
        return pieces.begin();
    }

    inline typename Pieces::iterator end()
    {
        return pieces.end();
    }

    inline void reserve(const int &n)
    {
        pieces.reserve(n);
        return;
    }

    inline void emplace_back(const Piece<D> &piece)
    {
        pieces.emplace_back(piece);
        return;
    }

    inline void emplace_back(const double &dur,
                             const typename Piece<D>::CoefficientMat &cMat)
    {
        pieces.emplace_back(dur, cMat);
        return;
    }

    inline void append(const Trajectory<D> &traj)
    {
        pieces.insert(pieces.end(), traj.begin(), traj.end());
        return;
    }

    inline int locatePieceIdx(double &t) const
    {
        int N = getPieceNum();
        int idx;
        double dur;
        for (idx = 0;
             idx < N &&
             t > (dur = pieces[idx].getDuration());
             idx++)
        {
            t -= dur;
        }
        if (idx == N)
        {
            idx--;
            t += pieces[idx].getDuration();
        }
        return idx;
    }

    inline Eigen::Vector3d getPos(double t) const
    {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getPos(t);
    }

    inline Eigen::Vector3d getVel(double t) const
    {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getVel(t);
    }

    inline Eigen::Vector3d getAcc(double t) const
    {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getAcc(t);
    }

    inline Eigen::Vector3d getJer(double t) const
    {
        int pieceIdx = locatePieceIdx(t);
        return pieces[pieceIdx].getJer(t);
    }

    inline Eigen::Vector3d getJuncPos(int juncIdx) const
    {
        if (juncIdx != getPieceNum())
        {
            return pieces[juncIdx].getCoeffMat().col(D);
        }
        else
        {
            return pieces[juncIdx - 1].getPos(pieces[juncIdx - 1].getDuration());
        }
    }

    inline Eigen::Vector3d getJuncVel(int juncIdx) const
    {
        if (juncIdx != getPieceNum())
        {
            return pieces[juncIdx].getCoeffMat().col(D - 1);
        }
        else
        {
            return pieces[juncIdx - 1].getVel(pieces[juncIdx - 1].getDuration());
        }
    }

    inline Eigen::Vector3d getJuncAcc(int juncIdx) const
    {
        if (juncIdx != getPieceNum())
        {
            return pieces[juncIdx].getCoeffMat().col(D - 2) * 2.0;
        }
        else
        {
            return pieces[juncIdx - 1].getAcc(pieces[juncIdx - 1].getDuration());
        }
    }
};

#endif