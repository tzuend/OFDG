// New experiment; include the preserved baseline without altering its source.
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreturn-type"
#define main preserved_direct_study_main
#include "direct_derivatives.cpp"
#undef main
#pragma clang diagnostic pop
#include "fe_mapping_jets.hpp"
struct FEElementData : direct::ElementData {std::unique_ptr<fe_mapping::GeometryField> field;};
int Run2D(int argc,char **argv) {

    if(argc!=8){std::cerr<<"family n p amplitude basis map extra\n";return 2;}
    std::string family=argv[1],basis=argv[5],mapping=argv[6];int n=std::stoi(argv[2]),p=std::stoi(argv[3]),extra=std::stoi(argv[7]);double amplitude=std::stod(argv[4]),h=1./n;
    bool triangle=family=="tri";
    Mesh mesh=Mesh::MakeCartesian2D(n,n,triangle?Element::TRIANGLE:Element::QUADRILATERAL,true);
    std::vector<FEElementData> data(mesh.GetNE());
    double initialization_seconds=0,evaluation_seconds=0,geometry_derivative_residual=0,geometry_jacobian_residual=0,fe_inverse_residual=0,method_agreement=0,fe_first_residual=0,face_geometry_residual=0;
    for(int e=0;e<mesh.GetNE();++e){auto *tr=mesh.GetElementTransformation(e);IntegrationPoint ip;ip.Set2(0,0);Vector x;tr->Transform(ip,x);tr->SetIntPoint(&ip);auto &j=tr->Jacobian();data[e].geometry={x(0),x(1),j(0,0),j(0,1),j(1,0),j(1,1),amplitude,mapping=="coupled"};}
    mesh.SetCurvature(4);mesh.Transform([&](const Vector&s,Vector&x){x=s;if(mapping=="coupled"){double b=amplitude*s(0)*(1-s(0))*s(1)*(1-s(1));x(0)+=b;x(1)+=.7*b;}else{x(0)+=amplitude*s(0)*(1-s(0));x(1)+=.7*amplitude*s(1)*(1-s(1));}});
    auto *first_iso=dynamic_cast<IsoparametricTransformation*>(mesh.GetElementTransformation(0));
    MFEM_VERIFY(first_iso,"Polynomial mesh geometry required");
    double setup_start=fe_mapping::Seconds();
    auto geometry_basis=std::make_shared<fe_mapping::PolynomialBasis>(*first_iso->GetFE());
    for(int e=0;e<mesh.GetNE();++e){auto *iso=dynamic_cast<IsoparametricTransformation*>(mesh.GetElementTransformation(e));
        MFEM_VERIFY(iso && iso->GetFE()==first_iso->GetFE(),"Uniform geometry basis required in this fixture");
        data[e].field=std::make_unique<fe_mapping::GeometryField>(*iso,geometry_basis);
    }
    initialization_seconds+=fe_mapping::Seconds()-setup_start;
    int bt=basis=="gll"?BasisType::GaussLobatto:(basis=="gl"?BasisType::GaussLegendre:BasisType::Positive);
    DG_FECollection fec(p,2,bt);FiniteElementSpace space(&mesh,&fec);direct::PolynomialBasis poly(*space.GetFE(0),triangle);
    std::vector<int> fields=mapping=="coupled"?std::vector<int>{0,1,2,3}:std::vector<int>{4,3};if(amplitude==0)fields.push_back(5);
    std::vector<std::string> words={"x","y","xx","xy","yy"};if(p>=3)for(auto s:{"xxx","xxy","xyy","yyy"})words.push_back(s);
    std::vector<direct::Row> rows;
    for(int f:fields)for(auto word:words)for(int mode=0;mode<3;++mode){int nx=std::count(word.begin(),word.end(),'x');rows.push_back({f,nx,int(word.size())-nx,mode,word});}
    double input[6]={},inverse_residual=0,representation_residual=0,first_derivative_residual=0,geometry_residual=0,minjac=1e100,face_measure=0,interior_measure=0;int qb=0;
    auto evaluate=[&](int e,const IntegrationPoint &ip,Vector &shape,std::array<direct::Jet,6>&values,std::array<direct::Jet,6>&fe_values) {
        double eval_start=fe_mapping::Seconds();
        const auto &ed=data[e];auto *fe=space.GetFE(e);auto *map=space.GetElementTransformation(e);fe->CalcShape(ip,shape);map->SetIntPoint(&ip);
        auto fe_inverse=fe_mapping::Inverse([&](const fe_mapping::Jets &r){return ed.field->Map(r);},ip,2,h,fe_inverse_residual);
        direct::Jet fr,fs;
        for(int i=0;i<direct::count;++i){int k=fe_mapping::Index(direct::ax[i],direct::ay[i]);fr.c[i]=fe_inverse[0].c[k];fs.c[i]=fe_inverse[1].c[k];}
        auto fe_monomials=poly.Evaluate(fr,fs,p);
        auto fe_geometry=ed.field->Map(fe_mapping::Coordinates(ip,2,true));
        direct::Jet ar(ip.x),as(ip.y);ar.c[1]=1;as.c[2]=1;auto analytic_geometry=ed.geometry.Map(ar,as);
        for(int d=0;d<2;++d)for(int i=0;i<direct::count;++i){int k=fe_mapping::Index(direct::ax[i],direct::ay[i]);
            geometry_derivative_residual=std::max(geometry_derivative_residual,std::abs(fe_geometry[d].c[k]-analytic_geometry[d].c[i])/h);}
        for(int d=0;d<2;++d)for(int k=0;k<2;++k)geometry_jacobian_residual=std::max(geometry_jacobian_residual,std::abs(fe_geometry[d].c[k+1]-map->Jacobian()(d,k))/h);
        auto rs=ed.geometry.InverseJet(ip.x,ip.y,h,inverse_residual);auto monomials=poly.Evaluate(rs[0],rs[1],p);
        DenseMatrix ds(fe->GetDof(),2);fe->CalcPhysDShape(*map,ds);Vector x;map->Transform(ip,x);
        auto xy=ed.geometry.Map(direct::Jet(ip.x),direct::Jet(ip.y));
        geometry_residual=std::max(geometry_residual,std::max(std::abs(xy[0].c[0]-x(0)),std::abs(xy[1].c[0]-x(1)))/h);
        for(int f:fields){fe_values[f]=direct::Jet();for(int j=0;j<int(fe_monomials.size());++j)fe_values[f]=fe_values[f]+ed.polynomial[f](j)*fe_monomials[j];
            values[f]=direct::Jet();for(int j=0;j<int(monomials.size());++j)values[f]=values[f]+ed.polynomial[f](j)*monomials[j];
            representation_residual=std::max(representation_residual,std::abs(values[f].c[0]-shape*ed.original[f]));
            for(int i=1;i<direct::count;++i){double fct=direct::Factorial(direct::ax[i])*direct::Factorial(direct::ay[i]);method_agreement=std::max(method_agreement,fct*std::abs(fe_values[f].c[i]-values[f].c[i])/std::max(1.,std::abs(fct*values[f].c[i])));}
            for(int axis=0;axis<2;++axis){double v=0;for(int j=0;j<fe->GetDof();++j)v+=ds(j,axis)*ed.original[f](j);fe_first_residual=std::max(fe_first_residual,std::abs(fe_values[f].c[axis+1]/h-v)/std::max(1.,std::abs(v)));first_derivative_residual=std::max(first_derivative_residual,std::abs(values[f].c[axis+1]/h-v)/std::max(1.,std::abs(v)));}
        }
        evaluation_seconds+=fe_mapping::Seconds()-eval_start;
    };
    auto derivative=[&](const direct::Jet &value,const direct::Row &r){return value.c[direct::Index(r.nx,r.ny)]*direct::Factorial(r.nx)*direct::Factorial(r.ny)/std::pow(h,r.nx+r.ny);};
    for(int e=0;e<space.GetNE();++e){double init_start=fe_mapping::Seconds();auto *fe=space.GetFE(e);auto *map=space.GetElementTransformation(e);int nd=fe->GetDof();auto &ed=data[e];qb=ofdg::detail::CurvedQuadratureOrder(p,*map);
        // Exactly the existing input projection and production derivative matrices.
        auto D=ofdg::detail::ProjectedPhysicalDerivatives(*fe,*map,qb);const auto &init=ofdg::detail::CurvedRule(fe->GetGeomType(),qb+16);DenseMatrix mass(nd);mass=0.;std::vector<Vector> rhs(6,Vector(nd));for(auto &v:rhs)v=0.;Vector shape(nd),x;
        for(int q=0;q<init.GetNPoints();++q){auto &ip=init.IntPoint(q);map->SetIntPoint(&ip);map->Transform(ip,x);fe->CalcShape(ip,shape);double w=ip.weight*map->Weight();AddMult_a_VVt(w,shape,mass);for(int f:fields)rhs[f].Add(w*Exact(f,x,0,0,amplitude),shape);}
        DenseMatrixInverse solve(mass);ed.original.resize(6,Vector(nd));ed.polynomial.resize(6,Vector(nd));ed.recursive.resize(rows.size());
        for(int f:fields){solve.Mult(rhs[f],ed.original[f]);poly.coefficients.Mult(ed.original[f],ed.polynomial[f]);}
        for(size_t i=0;i<rows.size();++i)if(rows[i].mode==0)ed.recursive[i]=Apply(D,ed.original[rows[i].field],rows[i].word);
        initialization_seconds+=fe_mapping::Seconds()-init_start;
        const auto &ir=ofdg::detail::CurvedRule(fe->GetGeomType(),qb+16+extra);
        for(int q=0;q<ir.GetNPoints();++q){auto &ip=ir.IntPoint(q);std::array<direct::Jet,6> values,fe_values;evaluate(e,ip,shape,values,fe_values);map->Transform(ip,x);double w=ip.weight*map->Weight();minjac=std::min(minjac,map->Jacobian().Det());
            for(int f:fields){double d=shape*ed.original[f]-Exact(f,x,0,0,amplitude);input[f]+=w*d*d;}
            for(size_t i=0;i<rows.size();++i){auto &r=rows[i];double exact=Exact(r.field,x,r.nx,r.ny,amplitude),dv=derivative(values[r.field],r),v=r.mode==2?derivative(fe_values[r.field],r):(r.mode==1?dv:shape*ed.recursive[i]),err=v-exact;r.error+=w*err*err;r.truth+=w*exact*exact;r.linf=std::max(r.linf,std::abs(err));r.operatorerror+=w*(v-dv)*(v-dv);r.directnorm+=w*dv*dv;}
        }
    }
    for(int face=0;face<mesh.GetNumFaces();++face){auto *tr=mesh.GetFaceElementTransformations(face);if(!tr)continue;const auto &ir=IntRules.Get(tr->GetGeometryType(),qb+16+extra);
        for(int q=0;q<ir.GetNPoints();++q){auto &ip=ir.IntPoint(q);tr->SetAllIntPoints(&ip);Vector x;tr->Face->Transform(ip,x);double w=ip.weight*tr->Face->Weight();std::vector<double> left(rows.size());
            for(int side=0;side<2;++side){int e=side?tr->Elem2No:tr->Elem1No;if(e<0)continue;face_measure+=w;const auto &local=side?tr->GetElement2IntPoint():tr->GetElement1IntPoint();Vector shape(space.GetFE(e)->GetDof());std::array<direct::Jet,6> values,fe_values;evaluate(e,local,shape,values,fe_values);
                auto fx=data[e].field->Map(fe_mapping::Coordinates(local,2));for(int d=0;d<2;++d)face_geometry_residual=std::max(face_geometry_residual,std::abs(fx[d].c[0]-x(d))/h);
                for(size_t i=0;i<rows.size();++i){auto &r=rows[i];double v=r.mode==2?derivative(fe_values[r.field],r):(r.mode==1?derivative(values[r.field],r):shape*data[e].recursive[i]),t=Exact(r.field,x,r.nx,r.ny,amplitude);r.face+=w*(v-t)*(v-t);r.facetruth+=w*t*t;if(side){r.jump+=w*(v-left[i])*(v-left[i]);r.jumpabs+=w*std::abs(v-left[i]);}else left[i]=v;}
                if(side)interior_measure+=w;
            }
        }
    }
    if(geometry_derivative_residual>1e-10||geometry_jacobian_residual>1e-10||fe_inverse_residual>1e-10||method_agreement>1e-8||fe_first_residual>1e-8||face_geometry_residual>1e-10){std::cerr<<"FE controls failed "<<geometry_derivative_residual<<' '<<geometry_jacobian_residual<<' '<<fe_inverse_residual<<' '<<method_agreement<<' '<<fe_first_residual<<' '<<face_geometry_residual<<'\n';return 4;}
    if(minjac<=0||inverse_residual>1e-9||representation_residual>1e-8||first_derivative_residual>1e-7||geometry_residual>1e-10){std::cerr<<"Control failed: "<<minjac<<' '<<inverse_residual<<' '<<representation_residual<<' '<<first_derivative_residual<<' '<<geometry_residual<<'\n';return 3;}
    std::cout<<"family,n,p,amplitude,basis,map,extra,field,derivative,method,error_l2,relative_l2,sampled_linf,truth_l2,face_rms,face_relative_l2,jump_rms,jump_mean_abs,input_error_l2,operator_error_l2,operator_relative_l2,inverse_residual,representation_residual,first_derivative_residual,geometry_residual,min_jacobian,elements,geometry_derivative_residual,geometry_jacobian_residual,fe_inverse_residual,method_agreement,fe_first_residual,face_geometry_residual,initialization_seconds,evaluation_seconds\n"<<std::setprecision(17);
    const char *names[]={"sine","shortwave","exponential","constant","represented_quadratic","physical_quadratic"};
    for(const auto &r:rows)std::cout<<family<<','<<n<<','<<p<<','<<amplitude<<','<<basis<<','<<mapping<<','<<extra<<','<<names[r.field]<<','<<r.word<<','<<(r.mode==2?"fe_geometry":(r.mode==1?"direct":"recursive"))<<','<<std::sqrt(r.error)<<','<<(r.truth>1e-24?std::sqrt(r.error/r.truth):NAN)<<','<<r.linf<<','<<std::sqrt(r.truth)<<','<<std::sqrt(r.face/face_measure)<<','<<(r.facetruth>1e-24?std::sqrt(r.face/r.facetruth):NAN)<<','<<std::sqrt(r.jump/interior_measure)<<','<<r.jumpabs/interior_measure<<','<<std::sqrt(input[r.field])<<','<<std::sqrt(r.operatorerror)<<','<<(r.directnorm>1e-24?std::sqrt(r.operatorerror/r.directnorm):NAN)<<','<<inverse_residual<<','<<representation_residual<<','<<first_derivative_residual<<','<<geometry_residual<<','<<minjac<<','<<space.GetNE()<<','<<geometry_derivative_residual<<','<<geometry_jacobian_residual<<','<<fe_inverse_residual<<','<<method_agreement<<','<<fe_first_residual<<','<<face_geometry_residual<<','<<initialization_seconds<<','<<evaluation_seconds<<'\n';
    return 0;
}

int main(int argc,char **argv){Mpi::Init(argc,argv);return Run2D(argc,argv);}
