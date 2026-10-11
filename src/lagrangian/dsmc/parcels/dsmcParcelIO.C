/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "dsmcParcel.H"
#include "IOstreams.H"
#include "IOField.H"
#include "dsmcCloud.H"
#include "PstreamReduceOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dsmcParcel::dsmcParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields
)
:
    particle(mesh, is, readFields),
    U_(vector::zero),
    RWF_(1.0),
    ERot_(0.0),
    ELevel_(0),
    typeId_(-1),
    newParcel_(-1),
    //tracked_(nullptr), // TODO TO BE REINTRODUCED
    //tracked_(0,0.0,vector::zero),
    classification_(0),
    stuck_(nullptr),
    vibLevel_(0)
{
    //dsmcParcel::TrackedParcel tP = dsmcParcel::TrackedParcel(); // TODO uncomment vincent  11/05/2018
    dsmcParcel::StuckParcel sP = dsmcParcel::StuckParcel();

    if (readFields)
    {
        if (is.format() == IOstream::ASCII)
        {
            is >> U_;
            RWF_ = readScalar(is);
            ERot_ = readScalar(is);
            ELevel_ = readLabel(is);
            typeId_ = readLabel(is);
            newParcel_ = readLabel(is);
            is >> tracked_; //tP; // TODO VINCENT  11/05/2018
            classification_ = readLabel(is);
            is >> sP;
            is >> vibLevel_;
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&U_),
                sizeof(U_)
                + sizeof(RWF_)
                + sizeof(ERot_)
                + sizeof(ELevel_)
                + sizeof(typeId_)
                + sizeof(newParcel_)
                + sizeof(classification_)
            );
            is >> tracked_; //tP;  // TODO VINCENT  11/05/2018
            is >> sP;
            is >> vibLevel_;
        }
    }

    if (false) // tP.tracked())
    {
        /*tracked_ = new dsmcParcel::TrackedParcel
            (
                tP.tracked(),
                tP.inPatchId(),
                tP.storePositions(),
                tP.initialTime(),
                tP.initialPosition(),
                tP.distanceTravelledVector()//,
                //tP.parcelTrajectory()
            );*/ // TODO TO BE REINTRODUCED

        /*tracked_ = dsmcParcel::TrackedParcel
            (
                tP.tracked(),
                tP.inPatchId(),
                tP.storePositions(),
                tP.initialTime(),
                tP.initialPosition(),
                tP.distanceTravelledVector()//,
                //tP.parcelTrajectory()
            );  */ // TODO uncomment vincent  11/05/2018
    }

    if (sP.wallTemperature()[0] != 0.0)
    {
        stuck_ = new dsmcParcel::StuckParcel
            (
                sP.wallTemperature(),
                sP.wallVectors()
            );
    }

    // Check state of Istream
    is.check
    (
        "Foam::dsmcParcel::dsmcParcel"
        "(const Cloud<dsmcParcel>& cloud, Foam::Istream&), bool"
    );
}


void Foam::dsmcParcel::readFields(Cloud<dsmcParcel>& c)
{
    if (!c.size())
    {
        return;
    }

    particle::readFields(c);

    IOField<vector> U(c.fieldIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    IOField<scalar> RWF
    (
        c.fieldIOobject
        (
            "radialWeight",
            IOobject::READ_IF_PRESENT
        ),
        scalarField(c.size(), 1.0)
    );
    c.checkFieldIOobject(c, RWF);

    IOField<scalar> ERot
    (
        c.fieldIOobject
        (
            "ERot",
            IOobject::READ_IF_PRESENT
        ),
        scalarField(c.size(), 0.0)
    );
    c.checkFieldIOobject(c, ERot);

    IOField<label> ELevel
    (
        c.fieldIOobject
        (
            "ELevel",
            IOobject::READ_IF_PRESENT
        ),
        labelField(c.size(), 0)
    );
    c.checkFieldIOobject(c, ELevel);

    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::MUST_READ));
    c.checkFieldIOobject(c, typeId);

    IOField<label> newParcel(c.fieldIOobject("newParcel", IOobject::MUST_READ));
    c.checkFieldIOobject(c, newParcel);

    IOField<label> classification(c.fieldIOobject("classification", IOobject::MUST_READ));
    c.checkFieldIOobject(c, classification);

    IOField<label> stuckToWall
    (
        c.fieldIOobject
        (
            "stuckToWall",
            IOobject::READ_IF_PRESENT
        ),
        labelField(c.size(), 0)
    );
    c.checkFieldIOobject(c, stuckToWall);

    IOField<scalarField> wallTemperature
    (
        c.fieldIOobject
        (
            "wallTemperature",
            IOobject::READ_IF_PRESENT
        )
    );

    if (wallTemperature.size() != c.size())
    {
        wallTemperature.setSize(c.size());
        forAll(wallTemperature, i)
        {
            wallTemperature[i] = scalarField(4, 0.0);
        }
    }

    c.checkFieldIOobject(c, wallTemperature);

    IOField<vectorField> wallVectors
    (
        c.fieldIOobject
        (
            "wallVectors",
            IOobject::READ_IF_PRESENT
        )
    );

    if (wallVectors.size() != c.size())
    {
        wallVectors.setSize(c.size());
        forAll(wallVectors, i)
        {
            wallVectors[i] = vectorField(4, vector::zero);
        }
    }

    c.checkFieldIOobject(c, wallVectors);

    IOField<label> isTracked
    (
        c.fieldIOobject
        (
            "isTracked",
            IOobject::READ_IF_PRESENT
        ),
        labelField(c.size(), 0)
    );
    c.checkFieldIOobject(c, isTracked);

    IOField<label> inPatchId
    (
        c.fieldIOobject
        (
            "inPatchId",
            IOobject::READ_IF_PRESENT
        ),
        labelField(c.size(), -1)
    );
    c.checkFieldIOobject(c, inPatchId);

    IOField<scalar> tracerInitialTime
    (
        c.fieldIOobject
        (
            "tracerInitialTime",
            IOobject::READ_IF_PRESENT
        ),
        scalarField(c.size(), 0.0)
    );
    c.checkFieldIOobject(c, tracerInitialTime);

    IOField<vector> tracerInitialPosition
    (
        c.fieldIOobject
        (
            "tracerInitialPosition",
            IOobject::READ_IF_PRESENT
        ),
        vectorField(c.size(), vector::zero)
    );
    c.checkFieldIOobject(c, tracerInitialPosition);

    IOField<labelField> vibLevel
    (
        c.fieldIOobject
        (
            "vibLevel",
            IOobject::READ_IF_PRESENT
        )
    );

    if (vibLevel.size() != c.size())
    {
        vibLevel.setSize(c.size());
        forAll(vibLevel, i)
        {
            vibLevel[i].setSize(0);
        }
    }

    c.checkFieldIOobject(c, vibLevel);

    label i = 0;
    forAllIter(dsmcCloud, c, iter)
    {
        dsmcParcel& p = iter();

        p.U_ = U[i];
        p.RWF_ = RWF[i];
        p.ERot_ = ERot[i];
        p.ELevel_ = ELevel[i];
        p.typeId_ = typeId[i];
        p.newParcel_ = newParcel[i];
        p.classification_ = classification[i];

        if (stuckToWall[i])
        {
            p.setStuck(wallTemperature[i], wallVectors[i]);
        }

        if (isTracked[i])
        {
            p.setTracked
            (
                isTracked[i],
                inPatchId[i],
                tracerInitialTime[i],
                tracerInitialPosition[i]
            );
        }

        p.vibLevel_ = vibLevel[i];

        i++;
    }
}


void Foam::dsmcParcel::writeFields(const Cloud<dsmcParcel>& c)
{
    particle::writeFields(c);

    const label np = c.size();

    IOField<vector> U(c.fieldIOobject("U", IOobject::NO_READ), np);
    IOField<scalar> RWF(c.fieldIOobject("radialWeight", IOobject::NO_READ), np);
    IOField<scalar> ERot(c.fieldIOobject("ERot", IOobject::NO_READ), np);
    IOField<labelField> vibLevel(c.fieldIOobject("vibLevel", IOobject::NO_READ), np);
    IOField<label> ELevel(c.fieldIOobject("ELevel", IOobject::NO_READ), np);
    IOField<label> typeId(c.fieldIOobject("typeId", IOobject::NO_READ), np);
    IOField<label> newParcel(c.fieldIOobject("newParcel", IOobject::NO_READ), np);
    IOField<label> classification(c.fieldIOobject("classification", IOobject::NO_READ), np);

    IOField<label> stuckToWall(c.fieldIOobject("stuckToWall", IOobject::NO_READ), np);
    IOField<scalarField> wallTemperature(c.fieldIOobject("wallTemperature", IOobject::NO_READ), np);
    IOField<vectorField> wallVectors(c.fieldIOobject("wallVectors", IOobject::NO_READ), np);

    IOField<label> isTracked(c.fieldIOobject("isTracked", IOobject::NO_READ), np);
    IOField<label> inPatchId(c.fieldIOobject("inPatchId", IOobject::NO_READ), np);
    IOField<scalar> tracerInitialTime(c.fieldIOobject("tracerInitialTime", IOobject::NO_READ), np);
    IOField<vector> tracerInitialPosition(c.fieldIOobject("tracerInitialPosition", IOobject::NO_READ), np);
    IOField<vector> tracerCurrentPosition(c.fieldIOobject("tracerCurrentPosition", IOobject::NO_READ), np);
    IOField<vector> tracerDistanceTravelled(c.fieldIOobject("tracerDistanceTravelled", IOobject::NO_READ), np);

    label i = 0;
    forAllConstIter(dsmcCloud, c, iter)
    {
        const dsmcParcel& p = iter();

        U[i] = p.U();
        RWF[i] = p.RWF();
        ERot[i] = p.ERot();
        vibLevel[i] = p.vibLevel();
        ELevel[i] = p.ELevel();
        typeId[i] = p.typeId();
        newParcel[i] = p.newParcel();
        classification[i] = p.classification();

        stuckToWall[i] = p.isStuck();
        if (stuckToWall[i])
        {
            wallTemperature[i] = p.stuck().wallTemperature();
            wallVectors[i] = p.stuck().wallVectors();
        }

        isTracked[i] = p.isTracked();
        if (isTracked[i])
        {
            inPatchId[i] = p.tracked().inPatchId();
            tracerInitialTime[i] = p.tracked().initialTime();
            tracerInitialPosition[i] = p.tracked().initialPosition();
            tracerCurrentPosition[i] = p.tracked().currentPosition();
            tracerDistanceTravelled[i] = p.tracked().distanceTravelledVector();
        }
        else
        {
            inPatchId[i] = -1;
            tracerInitialTime[i] = 0;
            tracerInitialPosition[i] = vector::zero;
            tracerCurrentPosition[i] = vector::zero;
            tracerDistanceTravelled[i] = vector::zero;
        }

        i++;
    }

    U.write();

    // 这里只能使用本 rank 局部的 radialWeight 偏离量。
    // 原因是 OpenFOAM 的 lagrangian 写场路径在并行空云 rank 上不一定会
    // 调用到 dsmcParcel::writeFields()。如果某个 rank 有 parcel 并进入这里，
    // 另一个 rank 没有 parcel 且没有进入这里，那么任何全局 reduce 都会造成
    // MPI 等待不闭合，表现为求解器在写完 positions/U 后卡住。
    scalar maxLocalRWFDeviation = 0.0;

    forAll(RWF, i)
    {
        // 对当前 rank 上每个 parcel 计算
        //
        //     |RWF_i - 1|
        //
        // 只要本 rank 有任意 parcel 的统计权重偏离 1，就在本 rank 写出
        // radialWeight。其它空 rank 或全为 1 的 rank 可以不写该字段；
        // readFields() 对缺失 radialWeight 的默认值本来就是 1。
        maxLocalRWFDeviation =
            max(maxLocalRWFDeviation, mag(RWF[i] - scalar(1.0)));
    }

    if (maxLocalRWFDeviation > SMALL)
    {
        RWF.write();
    }

    scalar maxLocalERot = 0.0;

    forAll(ERot, i)
    {
        // 与 radialWeight 相同，这里也不能用 gMax(ERot)。
        // 并行写场时空云 rank 可能根本不会进入 writeFields()，全局归约会
        // 因参与 rank 不一致而死锁。局部判据的重启语义仍然正确：
        // 没写 ERot 的 rank 会在 readFields() 中使用默认 ERot=0。
        maxLocalERot = max(maxLocalERot, ERot[i]);
    }

    if (maxLocalERot > 0.0)
    {
        //- there is at least one molecule
        ERot.write();
    }

    label maxLocalELevel = 0;

    forAll(ELevel, i)
    {
        // 电子能级字段同样采用本 rank 局部判据。
        // 若本 rank 全部为 0，则缺失 ELevel 字段会在读入时回到默认 0。
        maxLocalELevel = max(maxLocalELevel, ELevel[i]);
    }

    if (maxLocalELevel > 0)
    {
        //- the electronic mode is activated
        ELevel.write();
    }

    typeId.write();
    newParcel.write();
    classification.write();

    label nLocalStuckToWall = 0;

    forAll(stuckToWall, i)
    {
        // stuckToWall 是可选字段，只有本 rank 真实存在粘附粒子时才需要写。
        // 不能跨 rank 求和，否则空云 rank 未参与时会卡死。
        nLocalStuckToWall += stuckToWall[i];
    }

    if (nLocalStuckToWall > 0)
    {
        //- there is at least one stickingWallPatch with a particle stuck on it
        stuckToWall.write();
        wallTemperature.write();
        wallVectors.write();
    }

    label nLocalTracked = 0;

    forAll(isTracked, i)
    {
        // tracer 字段也是可选字段，局部存在 tracked parcel 才写。
        // 其它 rank 缺失该字段时，readFields() 会给出 isTracked=0 等默认值。
        nLocalTracked += isTracked[i];
    }

    if (nLocalTracked > 0)
    {
        //- there is at least one tracked parcel in the domain
        isTracked.write();
        inPatchId.write();
        tracerInitialTime.write();
        tracerInitialPosition.write();
        tracerCurrentPosition.write();
        tracerDistanceTravelled.write();
    }

    vibLevel.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const dsmcParcel& p
)
{
    /*dsmcParcel::TrackedParcel tP = dsmcParcel::TrackedParcel();
    if (p.isTracked())
    {
        tP = p.tracked();
    }*/ //TODO vincent 11/05/2018

    dsmcParcel::StuckParcel sP = dsmcParcel::StuckParcel();
    if (p.isStuck())
    {
        sP = p.stuck();
    }

    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.U()
            << token::SPACE << p.RWF()
            << token::SPACE << p.ERot()
            << token::SPACE << p.ELevel()
            << token::SPACE << p.typeId()
            << token::SPACE << p.newParcel()
            << token::SPACE << p.tracked() //tP; // TODO vincent  11/05/2018
            << token::SPACE << p.classification()
            << token::SPACE << sP
            << token::SPACE << p.vibLevel();
    }
    else
    {
        os  << static_cast<const particle&>(p);

        os.write
        (
            reinterpret_cast<const char*>(&p.U_),
            sizeof(p.U())
            + sizeof(p.RWF())
            + sizeof(p.ERot())
            + sizeof(p.ELevel())
            + sizeof(p.typeId())
            + sizeof(p.newParcel())
            + sizeof(p.classification())
        );

        os << p.tracked(); //tP; // TODO vincent  11/05/2018
        os << sP;
        os << p.vibLevel();
    }

    // Check state of Ostream
    os.check
    (
        "Foam::Ostream& Foam::operator<<"
        "(Foam::Ostream&, const Foam::dsmcParcel&)"
    );

    return os;
}


// ************************************************************************* //
