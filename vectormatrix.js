class Vec{
    static transpose(m){
        let t=new Array();
        for(let i=0;i<m.length;i++){
            t[i]=new Array();
            for(let j=0;j<m[i].length;j++){
                t[i][j]=m[j][i];
            }
        }
        return t;
    }

    static vecMatrixMul(v,m){
        let t=new Array();
        t[0]=v[0]*m[0][0]+v[1]*m[0][1]+v[2]*m[0][2];
        t[1]=v[0]*m[1][0]+v[1]*m[1][1]+v[2]*m[1][2];
        t[2]=v[0]*m[2][0]+v[1]*m[2][1]+v[2]*m[2][2];

        return t;
    }

    //Subtracts two arrays (vectors), a-b
    static sub(a, b){
        let t = new Array();
        t[0] = a[0] - b[0];
        t[1] = a[1] - b[1];
        t[2] = a[2] - b[2];
        return t;
    }

    //Gets a rotation matrix about the x axis.  Angle R is in radians
    static getXRotationMatrix(r){
        let t=new Array();
        t[0]=new Array();
        t[1]=new Array();
        t[2]=new Array();

        t[0][0]=1;
        t[0][1]=0;
        t[0][2]=0;
        t[1][0]=0;
        t[1][1]=Math.cos(r);
        t[1][2]=Math.sin(r);
        t[2][0]=0;
        t[2][1]=-Math.sin(r);
        t[2][2]=Math.cos(r);

        return t;
    }

    //Gets a rotation matrix about the y axis.  Angle R is in radians
    static getYRotationMatrix(r){
        let t=new Array();
        t[0]=new Array();
        t[1]=new Array();
        t[2]=new Array();

        t[0][0]=Math.cos(r);
        t[0][1]=0;
        t[0][2]=-Math.sin(r);
        t[1][0]=0;
        t[1][1]=1;
        t[1][2]=0;
        t[2][0]=Math.sin(r);
        t[2][1]=0;
        t[2][2]=Math.cos(r);

        return t;
    }

    //Gets a rotation matrix about the z axis.  Angle R is in radians
    static getZRotationMatrix(r){
        let t=new Array();
        t[0]=new Array();
        t[1]=new Array();
        t[2]=new Array();

        t[0][0]=Math.cos(r);
        t[0][1]=Math.sin(r);
        t[0][2]=0;
        t[1][0]=-Math.sin(r);
        t[1][1]=Math.cos(r);
        t[1][2]=0;
        t[2][0]=0;
        t[2][1]=0;
        t[2][2]=1;

        return t;
    }

    //Matrix dot product
    static dot(a,b){
        let m=new Array();
        for(let i=0;i<a.length;i++){
            m[i]=new Array();
            for(let j=0;j<b[0].length;j++){
                let temp=0;
                for(let k=0;k<b.length;k++){
                    temp+=a[i][k]*b[k][j];
                }
                m[i][j]=temp;
            }
        }
        return m;
    }
}